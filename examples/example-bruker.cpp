#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <iterator>
#include <cstdio>
#include "fbi/tuple.h"
#include "fbi/tuplegenerator.h"
#include "fbi/fbi.h"
#include "fbi/connectedcomponents.h"
#include "centroid.h"
#include "splitter.h"
#include "xic.h"

#include "boost/date_time/posix_time/posix_time.hpp"


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
  typedef SetA<Centroid,1,0 > SetType;
  typedef SetType::ResultType ResultType;
  typedef unsigned int LabelType;

  ptime start = microsec_clock::universal_time();
  
  SNSplitter<SetType, LabelType> splitter(options);

  CentroidBoxGenerator gen(10,1);
  ResultType test = std::vector<std::vector<unsigned int> >();
  boost::function<ResultType(std::deque<Centroid>)> boostFunctor =
  boost::bind(
  &SetType::intersect<const std::deque<Centroid>&,const CentroidBoxGenerator&, const CentroidBoxGenerator&>
  , _1, gen, gen);
  
  std::cout << "Looking for overlaps" << std::endl;
  ResultType fullAdjList = splitter.findOverlaps(boostFunctor);
  ptime end = microsec_clock::universal_time();
  time_duration td = end - start;
  std::cout << "elapsed time in seconds: " 
    << td.total_seconds()
    << std::endl;
  std::cout << "finding connected components... ";
  std::vector<LabelType> labels;
  unsigned int nComponents = findConnectedComponents(fullAdjList, labels);
  fullAdjList.clear();
  ResultType().swap(fullAdjList);

  std::cout << "Rereading centroids file" << std::endl;
  std::vector<Centroid> centroids;
  std::ifstream ifs(options.inputfileName_);
  std::string str;
  while (std::getline(ifs, str)) {
    parseString(str, centroids, 0);
  }
  ifs.close();

  std::vector<unsigned int> counter(nComponents);
  std::cout << "Creating Xics" << std::endl;
  std::vector<Xic> xics = createXicVector(centroids.begin(), centroids.end(), labels.begin(), labels.end(), counter);

  std::ofstream ofs(options.outputfileName_.c_str());


  std::cout << "Writing xics to outputfile" << std::endl;
  for (std::vector<Xic>::size_type i = 0; i < xics.size(); ++i) {
    if (counter[i] >= options.minClusterSize_) {
      ofs << xics[i] << "\t" << counter[i] <<std::endl;
    }
  }
}

