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
  double mzWindowLow_;
  int snWindowLow_;
  double mzWindowHigh_;
  int snWindowHigh_;
  std::string inputfileName_;
  std::string outputfileName_;
};


struct Centroid 
{
  double mz_;
  double sn_;
  Centroid(const double& mz, const double & sn) 
    : mz_(mz), sn_(sn){}
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







int parseProgramOptions(int argc, char* argv[], ProgramOptions& options)
{
  namespace po = boost::program_options;
  std::string config_file;
  po::options_description generic("Generic options");
  generic.add_options()
    ("help", "Display this help message")
    ("config,c", po::value<std::string>(&config_file), "config file")
    ("inputfile,i", po::value<std::string>(&options.inputfileName_), "input file")
    ("outputfile,o", po::value<std::string>(&options.outputfileName_), "output file")
    ;

  po::options_description config("Allowed options");
  config.add_options()
      ("mzLow", po::value<double>(&options.mzWindowLow_)->default_value(
        std::numeric_limits<double>::min()),
        "Lower bound of the query window for m/z")
      ("snLow", po::value<int>(&options.snWindowLow_)->default_value(
        std::numeric_limits<int>::min()),
        "Lower bound of the query window for the scan number")
      ("mzHigh", po::value<double>(&options.mzWindowHigh_)->default_value(
        std::numeric_limits<double>::max()),
        "Upper bound of the query window for m/z")
      ("snHigh", po::value<int>(&options.snWindowHigh_)->default_value(
        std::numeric_limits<int>::max()),
         "Upper bound of the query window for the scan number");

  po::options_description cmdline_options("Options available via command line");
  cmdline_options.add(generic).add(config);

  po::options_description config_file_options(
    "Options available in the config file");
  config_file_options.add(config);

  po::options_description visible("Allowed options");
  visible.add(generic).add(config);
  
  po::positional_options_description p;
  p.add("inputfile", -1);
  
  po::variables_map vm;

  po::store(po::command_line_parser(argc, argv).options(
    cmdline_options).positional(p).run(), vm);
  po::notify(vm);    
  if (vm.count("config")){ 
    std::ifstream ifs(config_file.c_str());
    if (!ifs) {
      std::cout << "can't open config file: " << config_file << '\n';
      return 0;
    } else {
      std::cerr << "Using config file" << config_file << '\n';
      po::store(po::parse_config_file(ifs, config_file_options), vm);
      po::notify(vm);    
    }
  }

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
parseFile(ProgramOptions & options) {
  std::vector<Centroid> centroids;
  boost::iostreams::stream<boost::iostreams::mapped_file_source> ifs(options.inputfileName_);

  std::string str;
  double mz, rt, abundance, intensity;
  std::string pol, mode, mslevel, line, massrange;
  int unknown, numentries, sn;
  typedef boost::tokenizer<boost::char_separator<char> > Tokenizer;
  boost::char_separator<char> sep(",");
  sn = 1;

  while(getline(ifs, str)) {
    
    Tokenizer tokens(str, sep);
    Tokenizer::iterator it = tokens.begin();
    std::stringstream sstream;
    if (it == tokens.end()) continue;
    sstream << *(it++);
    sstream >> rt;
    sstream.clear();
    if (it == tokens.end()) continue;
    sstream << *(it++);
    sstream >> pol;
    sstream.clear();
    if (it == tokens.end()) continue;
    sstream << *(it++);
    sstream >> mode;
    sstream.clear();
    if (it == tokens.end()) continue;
    sstream << *(it++);
    sstream >> mslevel;
    sstream.clear();
    if (it == tokens.end()) continue;
    sstream << *(it++);
    sstream >> unknown;
    sstream.clear();
    if (it == tokens.end()) continue;
    sstream << *(it++);
    sstream >> line;
    sstream.clear();
    if (it == tokens.end()) continue;
    sstream << *(it++);
    sstream >> massrange;
    sstream.clear();
    if (it == tokens.end()) continue;
    sstream << *(it++);
    sstream >> numentries;
    sstream.clear();

    int numcheck = 0;
    while(it != tokens.end() && numcheck < numentries) {
      sstream << *(it++);
      sstream >> mz >> intensity;
      sstream.clear();
      numcheck++;
      centroids.push_back(Centroid(mz, sn));
    }
    ++sn;
} 
return centroids;

}





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
  std::vector<Centroid> centroids = parseFile(options);
  std::cout << centroids.size() << std::endl;

  ptime start = microsec_clock::universal_time();
  SetA<Centroid,1,0>::ResultType centroidResults = SetA<Centroid, 1,0>::intersect(centroids, CentroidBoxGenerator(10,0.51), CentroidBoxGenerator(10,0.51));
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
}
