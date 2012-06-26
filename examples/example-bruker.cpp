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
  
  ResultType fullAdjList = splitter.findOverlaps(boostFunctor);
  ptime end = microsec_clock::universal_time();
  time_duration td = end - start;
  std::cout << "elapsed time in seconds: " 
    << td.total_seconds()
    << std::endl;
  std::cout << "finding connected components... ";
  std::vector<LabelType> labels;
  unsigned int nComponents = findConnectedComponents(fullAdjList, labels);
  std::vector<unsigned int> counter(nComponents, 0);
  for (std::vector<LabelType>::size_type i = 0; i < labels.size(); ++i) {
    ++counter[labels[i] - 1];
  }
  unsigned int bigClusterCounter = 0;
  for (unsigned int i = 0; i < counter.size(); ++i) {
    if (counter[i] >= options.minClusterSize_) bigClusterCounter++;
  }

  std::cout << "Number of clusters after filtering: " << bigClusterCounter << std::endl;
  std::sort(counter.begin(), counter.end());
  std::reverse(counter.begin(), counter.end());
  
/*
  std::ofstream ofs((options.outputfileName_ + ".counts").c_str());
  ofs << "Count" << std::endl;
  typedef std::map<unsigned int, unsigned int>::const_iterator It; 
  for (It it = counter.begin(); it != counter.end(); ++it) {
    ofs << *it << std::endl;
  }
  std::ifstream ifs(options.inputfileName_.c_str());
  std::vector<Centroid> centroids;
  int sn = 0;
  std::string str;
  while (std::getline(ifs, str)) {
    parseString(str, centroids, sn);
  }


  typedef std::vector<LabelType>::const_iterator LabelIter;
  std::vector<Xic> xics;  
  xics.resize(bifgClusterCounter);
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
*/

}
