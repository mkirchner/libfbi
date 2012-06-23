#include <boost/tokenizer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <vector>
#include <fstream>
#include <cstdio>
#include "centroid.h"

struct ProgramOptions 
{
  std::string inputfileName_;
  std::string outputfileName_;
  unsigned int segmentsize_;
  unsigned int overlap_;
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



template <class SetType>
struct
SNSplitter {
  private:
    ProgramOptions options_;
    std::vector<unsigned int> breakpoints_;
  public:
    bool
      setOptions(const ProgramOptions & options) {
        return (options_ = options);
      }

    std::vector<Centroid>
      parseString(const std::string & str, unsigned int sn) {
        std::vector<Centroid> centroids;
        float mz, massrange_lo, massrange_hi, rt;
        char pol;
        char mode[101];
        char mslevel[101];
        char line[101];
        int unknown, numentries, intensity;
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
          ++sn;
        }
      }
  ResultType
  filterAdjList(const ResultType & adjList, unsigned int lowerLimit) {
    ResultType filteredAdjList;
    try {
      std::vector<LabelType> labels;
      LabelType nComponents = findConnectedComponents(filteredAdjList, labels); 
      std::vector<unsigned int> counter(*std::max_element(labels.begin(), labels.end()), 0);
      for (std::vector<unsigned int>::size_type i = 0; i < labels.size(); ++i) {
        ++counter[labels[i] - 1];
      }
      filteredAdjList.resize(adjList.size());
      for (typename ResultType::size_type i = 0; i < tempResult.size(); ++i) {
        if (counter[labels[i]-1] > rightoverlap_) { 
          std::copy(adjList[i].begin(), adjList[i].end(), std::back_inserter(filteredAdjList[i]));
          }
      }



    return filteredAdjList;
  }
    ResultType
    findOverlaps(boost::function<ResultType(std::vector<Centroid>)> intersectFunctor) {
       
        boost::iostreams::stream<boost::iostreams::file_source> ifs(options.inputfileName_);
        std::string str;
        std::deque<Centroid> centroids;
        unsigned int segmentCounter = 0;
        unsigned int centroidCounter = 0;
        ResultType fullAdjList;
        while(std::getline(ifs, str)) {
          unsigned int breakpoint = parseString(str, centroids, counter);
          segmentCounter++;
          if (segmentCounter % options.segmentsize_ == 0) {
            centroidCounter = centroids.size();
          }
          if (segmentCounter % options.segmentsize_ == overlap && segmentCounter > options.segmentsize_) {
            ResultType shortAdjList = intersectFunctor(centroids);
            ResultType filteredAdjList = filterAdjList(shortAdjList, labels);
            joinAdjLists(filteredAdjList, fullAdjList, centroidCounter);
          }
        }
        return fullAdjList;


    }
    
};
