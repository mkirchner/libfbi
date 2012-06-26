#include <boost/tokenizer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <vector>
#include <deque>
#include <fstream>
#include <cstdio>
#include "centroid.h"
#ifndef __LIBFBI_EXAMPLES_SPLITTER_H__
#define __LIBFBI_EXAMPLES_SPLITTER_H__
struct ProgramOptions 
{
  std::string inputfileName_;
  std::string outputfileName_;
  unsigned int segmentSize_;
  unsigned int minClusterSize_;
  unsigned int snWindowSize_;
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
      ("segmentSize", po::value<unsigned int>(&options.segmentSize_)->default_value(
        100),
        "The number of lines read before calculating adjList for each partition")
      ("minClusterSize", po::value<unsigned int>(&options.minClusterSize_)->default_value(
        5),
         "The size a cluster has to have at least to not be purged")
      ("snWindowSize", po::value<unsigned int>(&options.snWindowSize_)->default_value(1),
      "The maximum distance between two points in the same cluster in SN direction");


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
template <typename ContainerType>
unsigned int 
parseString(const std::string & str, ContainerType & centroids, int sn, bool onlyheader = false) {
  float mz=0, massrange_lo=0, massrange_hi=0, rt=0;
  char pol=0;
  char mode[101];
  char mslevel[101];
  char line[101];
  int unknown=0, numentries=0, intensity=0;
  typedef boost::tokenizer<boost::char_separator<char> > Tokenizer;
  boost::char_separator<char> sep(",");

  Tokenizer tokens(str, sep);
  Tokenizer::iterator it = tokens.begin();
  if (sscanf(str.c_str(), "%f,%c,%100[^,],%100[^,],%d,%100[^,],%f-%f,%u,", &rt, &pol, mode, mslevel, &unknown, line, &massrange_lo, &massrange_hi, &numentries) != 9 || onlyheader) {
    return numentries;
  }
  for (int i = 0; i < 8; ++i) ++it;
  while (it != tokens.end()) {
    std::string mz_int_pair(*(it++));
    if (sscanf(mz_int_pair.c_str(), "%f %u", &mz, &intensity) == 2) {
      centroids.push_back(Centroid(mz, sn, rt));
    }
  }
  return numentries;
}

template <class SetType, class LabelType = unsigned int>
struct
SNSplitter {
  private:
    ProgramOptions options_;
    unsigned int overlap_;
    unsigned int adjListCounter_;

  public:
    typedef typename SetType::ResultType ResultType;
    typedef typename std::vector<unsigned int>::size_type size_type;
    SNSplitter(const ProgramOptions& options) : options_(options) {

      overlap_ = std::max(options_.minClusterSize_ * options_.snWindowSize_ - 1, options_.snWindowSize_);
      adjListCounter_ = 0;
    }

    ResultType
      filterAdjList(ResultType & adjList) {
        std::vector<LabelType> labels;
        LabelType nComponents = findConnectedComponents(adjList, labels); 
        std::vector<unsigned int> counter(nComponents, 0);
        for (size_type i = 0; i < labels.size(); ++i) {
          ++counter[labels[i] - 1];
        }
        typedef typename ResultType::value_type InnerType;
        for (size_type i = 0; i < labels.size(); ++i) {
          if (counter[labels[i]-1] < options_.minClusterSize_) { 
            InnerType().swap(adjList[i]);
            //std::copy(adjList[i].begin(), adjList[i].end(), std::back_inserter(filteredAdjList[i]));
          }
        }
        return adjList;
      }
    ResultType
      joinAdjLists(const ResultType & filteredAdjList, ResultType & fullAdjList, size_type offset) {
        //fullAdjList.resize(filteredAdjList.size() + offset);
        size_type index;
        typedef typename ResultType::value_type InnerType;
        typename ResultType::value_type::const_iterator it1;
        for (index = 0; index < filteredAdjList.size(); ++index) {
          for (it1 = filteredAdjList[index].begin(); it1 != filteredAdjList[index].end(); ++it1) {
            fullAdjList[index + offset].insert(fullAdjList[index + offset].end(), typename ResultType::value_type::value_type(*it1 + offset));
            adjListCounter_++;
          }
          //InnerType().swap(filteredAdjList[index]);
        }

        return fullAdjList;
      }
    ResultType
      createShortAdjList(boost::function<ResultType(std::deque<Centroid>)> & intersectFunctor, const std::deque<Centroid> & centroids) {
        ResultType shortAdjList = intersectFunctor(centroids);

        filterAdjList(shortAdjList);

        return shortAdjList;
      }
    template <typename IteratorType>
      IteratorType
      makeUnique(IteratorType begin, IteratorType end) {
        for (;begin != end; ++begin) {
          typename ResultType::value_type & vec = *begin;
          std::sort(vec.begin(), vec.end());
          vec.resize(std::unique(vec.begin(), vec.end()) - vec.begin());
        }
        return begin;
      }


    ResultType
      findOverlaps(boost::function<ResultType(std::deque<Centroid>)> & intersectFunctor) {

        boost::iostreams::stream<boost::iostreams::file_source> ifs(options_.inputfileName_);
        std::string str;
        std::deque<Centroid> centroids;
        unsigned int segmentCounter = 0;
        unsigned int centroidCounter = 0;
        unsigned int oldCentroidCounter = 0;
        unsigned int nextCentroidCounter = 0;

        ResultType fullAdjList;
        unsigned int countAllCentroids = 0;
        while(std::getline(ifs, str)) {
          bool onlyheader = true;
          countAllCentroids += parseString(str, centroids, segmentCounter, onlyheader);
        }
        fullAdjList.resize(countAllCentroids);
        ifs.clear();
        ifs.seekg(0);
        while(std::getline(ifs, str)) {
          unsigned int numEntries = parseString(str, centroids, segmentCounter);
          segmentCounter++;
          centroidCounter += numEntries;
          if (segmentCounter % options_.segmentSize_ == 0) {
            nextCentroidCounter = centroidCounter;
          }

          if (segmentCounter % options_.segmentSize_ == overlap_ && segmentCounter >= options_.segmentSize_) {
            ResultType filteredAdjList = createShortAdjList(intersectFunctor, centroids);

            joinAdjLists(filteredAdjList, fullAdjList, oldCentroidCounter);
            makeUnique(fullAdjList.begin() + oldCentroidCounter, fullAdjList.begin() + centroidCounter);
            unsigned int numNewCentroids = nextCentroidCounter - oldCentroidCounter;
            oldCentroidCounter = nextCentroidCounter;
            centroids.erase(centroids.begin(), centroids.begin() + numNewCentroids);
            ResultType().swap(filteredAdjList); 

          }
        }
        ifs.close();
        ResultType filteredAdjList = createShortAdjList(intersectFunctor, centroids);
        joinAdjLists(filteredAdjList, fullAdjList, oldCentroidCounter);

        makeUnique(fullAdjList.begin()+oldCentroidCounter, fullAdjList.end());

        return fullAdjList;
      }

};
#endif
