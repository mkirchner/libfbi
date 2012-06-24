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



template <class SetType>
struct
SNSplitter {
  private:
    ProgramOptions options_;
    unsigned int overlap_;
    unsigned int adjListCounter_;
  public:
  typedef typename SetType::ResultType ResultType;
  typedef unsigned int LabelType;
    SNSplitter(const ProgramOptions& options) : options_(options) {
      overlap_ = options_.minClusterSize_ * options_.snWindowSize_ - 1;
      std::cout << overlap_ << "overlap";
      adjListCounter_ = 0;
    }
     unsigned int 
      parseString(const std::string & str, std::deque<Centroid> & centroids, unsigned int sn) {
        float mz=0, massrange_lo=0, massrange_hi=0, rt=0;
        char pol=0;
        char mode[101];
        char mslevel[101];
        char line[101];
        int unknown=0, numentries=0, intensity=0;
        typedef boost::tokenizer<boost::char_separator<char> > Tokenizer;
        boost::char_separator<char> sep(",");
        sn = 1;

        Tokenizer tokens(str, sep);
        Tokenizer::iterator it = tokens.begin();
        if (sscanf(str.c_str(), "%f,%c,%100[^,],%100[^,],%d,%100[^,],%f-%f,%u,", &rt, &pol, mode, mslevel, &unknown, line, &massrange_lo, &massrange_hi, &numentries) != 9) {
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
    ResultType
      filterAdjList(const ResultType & adjList) {
        ResultType filteredAdjList;
        std::vector<LabelType> labels;
        LabelType nComponents = findConnectedComponents(adjList, labels); 
        std::vector<unsigned int> counter(nComponents, 0);
        for (std::vector<LabelType>::size_type i = 0; i < labels.size(); ++i) {
          ++counter[labels[i] - 1];
        }
        //std::cout << "Counter" << counter;
        filteredAdjList.resize(adjList.size());
        for (std::vector<LabelType>::size_type i = 0; i < labels.size(); ++i) {
          if (counter[labels[i]-1] >= options_.minClusterSize_) { 
            std::copy(adjList[i].begin(), adjList[i].end(), std::back_inserter(filteredAdjList[i]));
          }
        }
        return filteredAdjList;
      }
      ResultType
      joinAdjLists(const ResultType & filteredAdjList, ResultType & fullAdjList, typename ResultType::size_type offset) {
        fullAdjList.resize(filteredAdjList.size() + offset);
        std::cout << "OFFSET" << offset << std::endl;
        typename ResultType::size_type index;
        typename ResultType::value_type::const_iterator it1;
        for (index = 0; index < filteredAdjList.size(); ++index) {
          for (it1 = filteredAdjList[index].begin(); it1 != filteredAdjList[index].end(); ++it1) {
            fullAdjList[index + offset].insert(fullAdjList[index + offset].end(), *it1 + offset);
          }
        }
        //std::cout << adjListCounter_ << "adjCounter" << std::endl;
        return fullAdjList;
      }
    ResultType
      createShortAdjList(boost::function<ResultType(std::deque<Centroid>)> & intersectFunctor, const std::deque<Centroid> & centroids) {
        ResultType shortAdjList = intersectFunctor(centroids);

        ResultType filteredAdjList = filterAdjList(shortAdjList);
        for (int i = 0; i < 10; ++i) {
        //std::cout << shortAdjList[i] << std::endl;
        }
        //std::cout << "FUBLA" << std::endl;
        for (int i = 0; i < 10; ++i) {
        //std::cout << filteredAdjList[i] << std::endl;
        }

        return filteredAdjList;
      }

    ResultType & makeUnique(ResultType & resultVector) {
      for (typename ResultType::size_type i = 0; i < resultVector.size(); ++i) {
        typename ResultType::value_type & vec = resultVector[i];
        std::sort(vec.begin(), vec.end());
        vec.resize(std::unique(vec.begin(), vec.end()) - vec.begin());
      }

      return resultVector;
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
        while(std::getline(ifs, str)) {
          unsigned int numEntries = parseString(str, centroids, segmentCounter);
          segmentCounter++;
          centroidCounter += numEntries;
          
          if (segmentCounter % options_.segmentSize_ == 0) {
            nextCentroidCounter = centroidCounter;
          }
          
          if (segmentCounter % options_.segmentSize_ == overlap_ && segmentCounter >= options_.segmentSize_) {
            std::cout << "pung" << std::endl;
            std::cout << "oldCentroidCounter" << oldCentroidCounter << std::endl;
            std::cout << "nextCentroidCounter" << nextCentroidCounter << std::endl;
            ResultType filteredAdjList = createShortAdjList(intersectFunctor, centroids);
            
            joinAdjLists(filteredAdjList, fullAdjList, oldCentroidCounter);
            std::cout << "centroid" << centroids.size() << std::endl;
            std::cout << "centroidCounter" << centroidCounter << std::endl;
            std::cout << "fullAdj" << fullAdjList.size() << std::endl;
            unsigned int numNewCentroids = nextCentroidCounter - oldCentroidCounter;
            oldCentroidCounter = nextCentroidCounter;
            centroids.erase(centroids.begin(), centroids.begin() + numNewCentroids);
            

            std::cout << "pong" << std::endl;
          }
          
        }
        
        ResultType filteredAdjList = createShortAdjList(intersectFunctor, centroids);
        joinAdjLists(filteredAdjList, fullAdjList, oldCentroidCounter);

        fullAdjList = makeUnique(fullAdjList);
        return fullAdjList;
      }

};
#endif
