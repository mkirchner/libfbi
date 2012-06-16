#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include "fbi/tuple.h"
#include "fbi/tuplegenerator.h"
#include "fbi/fbi.h"
#include "fbi/connectedcomponents.h"

#include "boost/date_time/posix_time/posix_time.hpp"

struct ProgramOptions 
{
  std::string inputfileName_;
  std::string outputfileName_;
  unsigned int segments_;
  unsigned int overlap_;
};


struct Centroid 
{
  double mz_;
  double sn_;
  double rt_;
  Centroid(const double& mz, const double & sn, const double & rt) 
    : mz_(mz), sn_(sn), rt_(rt){}
};



namespace fbi {

  template<>
  struct Traits<Centroid> : mpl::TraitsGenerator<float, float> {};

} //end namespace fbi

struct CentroidBoxGenerator
{
  template <size_t N>
    typename fbi::tuple_element<N, 
             typename fbi::Traits<Centroid>::key_type>::type 
               get(const Centroid&) const;

  double mzOffset_;
  double mzWindowPpm_;
  double rtOffset_;
  double rtWindow_;
  double snWindow_;

  CentroidBoxGenerator(double mzWindowPpm, double snWindow)
    : mzOffset_(0.0), mzWindowPpm_(mzWindowPpm), snWindow_(snWindow)
  {}
  CentroidBoxGenerator(double mzWindowPpm, double rtWindow, int snWindow)
    : mzOffset_(0.0), mzWindowPpm_(mzWindowPpm),
    rtOffset_(0.0), rtWindow_(rtWindow), snWindow_(snWindow)
  {}

  CentroidBoxGenerator(double mzOffset, double mzWindowPpm, 
      double rtOffset, double rtWindow, int snWindow)
    : mzOffset_(mzOffset), mzWindowPpm_(mzWindowPpm),
    rtOffset_(rtOffset), rtWindow_(rtWindow), snWindow_(snWindow)
  {}


};


template <>
std::pair<float, float>  
CentroidBoxGenerator::get<0>(const Centroid & centroid) const 
{
  return std::make_pair(
      mzOffset_ + centroid.mz_* (1 - mzWindowPpm_ * 1E-6), 
      mzOffset_ + centroid.mz_* (1 + mzWindowPpm_ * 1E-6));
}


template <>
std::pair<float, float>  
CentroidBoxGenerator::get<1>(const Centroid & centroid) const 
{
  return std::make_pair(
      centroid.sn_ - snWindow_, 
      centroid.sn_ + snWindow_);
}


struct Xic
{
  double mz_;
  double rt_;
  Xic(){}
  Xic(const double& mz, const double & rt) 
    : mz_(mz), rt_(rt){}
};




int parseProgramOptions(int argc, char* argv[], ProgramOptions& options)
{
  namespace po = boost::program_options;
  po::options_description generic("Generic options");
  generic.add_options()
    ("help", "Display this help message")
    ("inputfile,i", po::value<std::string>(&options.inputfileName_), "input file")
    ("outputfile,o", po::value<std::string>(&options.outputfileName_), "output file")
    ;

  po::options_description config("Allowed options");
  config.add_options()
      ("segments", po::value<unsigned int>(&options.segments_)->default_value(
        1),
        "Number of segments the data should be partitioned in before using libfbi")
      ("overlap", po::value<unsigned int>(&options.overlap_)->default_value(
        0),
         "Overlap in time-dimension taken into account to not have jumps");
  
  po::options_description cmdline_options("Options available via command line");
  cmdline_options.add(generic).add(config);

  po::options_description visible("Allowed options");
  visible.add(generic).add(config);

  po::positional_options_description p;
  p.add("inputfile", -1);

  po::variables_map vm;

  po::store(po::command_line_parser(argc, argv).options(
        cmdline_options).positional(p).run(), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    std::cout << visible << "\n";
    return 0;
  }

  if (!vm.count("inputfile")) {
    std::cerr << "InputFile needed" << '\n';
    std::cout << visible << '\n';
    return 0;
  }

  if (!vm.count("outputfile")) {
    options.outputfileName_ = options.inputfileName_ + std::string(".out");
  }
  return 1;
}


std::vector<Centroid>
parseFile(ProgramOptions & options, std::vector<std::vector<Centroid>::size_type> & breakpoints) {
  std::vector<Centroid> centroids;
  boost::iostreams::stream<boost::iostreams::mapped_file_source> ifs(options.inputfileName_);
  std::string str;
  float mz, massrange_lo, massrange_hi, rt;
  char pol;
  char mode[100];
  char mslevel[100];
  char line[100];
  int unknown, numentries, intensity, sn;
  typedef boost::tokenizer<boost::char_separator<char> > Tokenizer;
  boost::char_separator<char> sep(",");
  sn = 1;
  breakpoints.push_back(0);
  while(getline(ifs, str)) {

    Tokenizer tokens(str, sep);
    Tokenizer::iterator it = tokens.begin();
	if (sscanf(str.c_str(), "%f,%c,%100[^,],%100[^,],%d,%100[^,],%f-%f,%u,", &rt, &pol, &mode, &mslevel, &unknown, &line, &massrange_lo, &massrange_hi, &numentries) != 9) {
		continue;
	}
	for (int i = 0; i < 8; ++i) ++it;
	while (it != tokens.end()) {
		std::string mz_int_pair(*(it++));
		if (sscanf(mz_int_pair.c_str(), "%f %u", &mz, &intensity) == 2) {
			centroids.push_back(Centroid(mz, sn, rt));
		}
	}

    breakpoints.push_back(centroids.size());
    ++sn;
  } 
  return centroids;

}


//With the prior knowledge of having very small "time"-boxes, we can split up data beforehand,
//combining all of the results for the adjacency list afterwards.
template <typename ResultType>
struct SNSplitter{
  typedef typename ResultType::size_type size_type;

  //support for two generators, prototyping
  SNSplitter(const std::vector<Centroid>& centroids, const std::vector<std::vector<Centroid>::size_type> & breakpoints, CentroidBoxGenerator b1, CentroidBoxGenerator b2, unsigned int segments, unsigned int rightoverlap):
    centroids_(centroids), breakpoints_(breakpoints), segments_(segments), rightoverlap_(rightoverlap), b1_(b1), b2_(b2)
  {
    //if (segments_ == 0) segments_ = 1;
	  if (breakpoints.size() < segments_) segments_ = 1;
	
    offsets.resize(segments_);
    limits.resize(segments_);
    unsigned int stepsize = ((unsigned int)(breakpoints_.size()-1) / segments_);
    offsets[0] = 0;
    limits[0] = (unsigned int) breakpoints_[stepsize + rightoverlap];
    for (unsigned int i = 1; i < segments_ - 1; ++i) {
      offsets[i] = (unsigned int) breakpoints_[i*stepsize];
      limits[i] = (unsigned int) breakpoints_[(i+1)*stepsize+rightoverlap];
    }
	  offsets[segments_-1] = (unsigned int) breakpoints_[(segments_ - 1)*stepsize];
    limits[segments_-1] = (unsigned int)centroids_.size();    
  }


  ResultType & join(ResultType & completeResult, const ResultType & tempResult, unsigned int i) {
    completeResult.resize(limits[i]);
    typename ResultType::size_type index;
    typename ResultType::value_type::value_type currentoffset;
    typename ResultType::value_type::const_iterator it1;
    currentoffset = offsets[i];

    for (index = 0; index < tempResult.size(); ++index) {
      for (it1 = tempResult[index].begin(); it1 != tempResult[index].end(); ++it1) {
        completeResult[index + currentoffset].insert(completeResult[index + currentoffset].end(), *it1 + currentoffset);
      }
    }
    return completeResult;
  }

  ResultType & makeUnique(ResultType & resultVector) {
    for (typename ResultType::size_type i = 0; i < resultVector.size(); ++i) {
      typename ResultType::value_type & vec = resultVector[i];
      std::sort(vec.begin(), vec.end());
      vec.resize(std::unique(vec.begin(), vec.end()) - vec.begin());
    }

    return resultVector;
  }

  ResultType operator()() {
  using namespace fbi;
    ResultType completeResult;
	try {
    for (unsigned int i = 0; i < segments_; ++i) {
      std::vector<Centroid> shortList(centroids_.begin() + offsets[i], centroids_.begin() + limits[i]);
	  
      ResultType tempResult = SetA<Centroid, 1, 0>::intersect(shortList, b1_, b2_);
      completeResult = join(completeResult, tempResult, i);
    }
    completeResult = makeUnique(completeResult);
	}
	catch (const std::exception & e) {
		std::cout << e.what();
	}
    return completeResult;
  }


  const std::vector<Centroid> & centroids_;
  
  
  std::vector<size_type> breakpoints_; 

  unsigned int segments_;
  unsigned int rightoverlap_;
  CentroidBoxGenerator b1_, b2_;
  //These are the offsets to be added to the indices of later results, so the combination still makes sense.
  //Same as the first centroid belonging to a specific segment.
  std::vector<unsigned int> offsets;
  //The last centroid belonging to a specific segment.
  std::vector<unsigned int> limits;

};



int main(int argc, char * argv[]) {
  using namespace fbi;
  using namespace boost::posix_time;
#ifdef __LIBFBI_USE_MULTITHREADING__
  std::cout << "Multithreading enabled" << std::endl;
#endif
  ProgramOptions options;
  if (!parseProgramOptions(argc, argv, options)) {
    return 0;
  }
  std::vector<std::vector<Centroid>::size_type> breakpoints;
  std::vector<Centroid> centroids = parseFile(options, breakpoints);
  std::cout << centroids.size() << std::endl;

  ptime start = microsec_clock::universal_time();
  typedef SetA<Centroid,1,0 >::ResultType ResultType;
  //  SetA<Centroid,1,0>::ResultType centroidResults = SetA<Centroid, 1,0>::intersect(centroids, CentroidBoxGenerator(10,0.51), CentroidBoxGenerator(10,0.51));
  SNSplitter<ResultType> splitter(centroids,breakpoints, CentroidBoxGenerator(10,0.51), CentroidBoxGenerator(10,0.51), options.segments_,options.overlap_);
  ResultType centroidResults = splitter();


  ptime end = microsec_clock::universal_time();
  time_duration td = end - start;



  std::cout << "elapsed time in seconds: " 
    << td.total_seconds()
    << std::endl;

  std::cout << "finding connected components... ";
  typedef SetA<Centroid, 1,0>::IntType LabelType;
  std::vector<LabelType> labels;
  LabelType nComponents = findConnectedComponents(centroidResults, labels); 
  std::cout << nComponents << " components found." << std::endl;
  typedef std::vector<LabelType>::const_iterator LabelIter;
  std::vector<Xic> xics;  
  xics.resize(nComponents);
  std::vector<unsigned int> labelcounts;
  labelcounts.resize(nComponents, 0);
  for (unsigned int i = 0; i < labels.size(); ++i) {
	xics[labels[i] - 1].mz_ += centroids[i].mz_;
	xics[labels[i] - 1].rt_ += centroids[i].rt_;
	++labelcounts[labels[i] - 1];
  }
  for (unsigned int i = 0; i < xics.size(); ++i) {
	xics[i].mz_ /= labelcounts[i];
	xics[i].rt_ /= labelcounts[i];
  }

  std::sort(labelcounts.begin(), labelcounts.end());
  std::cout << "Biggest clusters: " << std::endl;
  for (unsigned int i = 1; i < 100; ++i) {
	  std::cout << i << ": " << labelcounts[labelcounts.size() - i] << std::endl;
  }
  std::ofstream ofs(options.outputfileName_.c_str());
  for (unsigned int i = 0; i < xics.size(); ++i) {
	  ofs << xics[i].mz_ << "\t" << xics[i].rt_ << std::endl;
  }
}
