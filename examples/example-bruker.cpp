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
#include <prettyprint.hpp>
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



  ptime start = microsec_clock::universal_time();
  
  SNSplitter<SetType> splitter(options);

  CentroidBoxGenerator gen(10,options.snWindowSize_);
  //ResultType (*intersectFunctor)
  //(std::vector<Centroid>, CentroidBoxGenerator, CentroidBoxGenerator); 
  
  //(ResultType
  //(SetType::*)
  //(const std::vector<Centroid>&,const CentroidBoxGenerator&, const CentroidBoxGenerator&) )
  //&SetType::intersect<const std::vector<Centroid>&,const CentroidBoxGenerator&, const CentroidBoxGenerator&>;
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
  std::vector<unsigned int> labels;
  unsigned int nComponents = findConnectedComponents(fullAdjList, labels);
  for (int i = 0; i < 100; ++i) {
  //std::cout << fullAdjList[i] << std::endl;
  }
  std::cout << nComponents << "nComponents" << std::endl;

/*

  std::cout << "finding connected components... ";
  typedef SetType::IntType LabelType;
  std::vector<LabelType> labels, filteredLabels;
  //LabelType nComponents = findConnectedComponents(centroidResults, labels); 
  //std::cout << nComponents << " components found." << std::endl;
  LabelType dComponents = findConnectedComponents(filteredResults, filteredLabels); 
  std::vector<unsigned int> labelCounts = countLabels(filteredLabels, dComponents);
  std::map<unsigned int, unsigned int> filteredCounts = 
    filterCounts(labelCounts, options.overlap_);

  std::cout << filteredCounts.size() << " components found." << std::endl;

  std::ofstream ofs(options.outputfileName_.c_str());
  ofs << "Label" << "\t" << "Count" << std::endl;
  typedef std::map<unsigned int, unsigned int>::const_iterator It; 
  for (It it = filteredCounts.begin(); it != filteredCounts.end(); ++it) {
    ofs << it->first << "\t" << it->second  << std::endl;
  }
*/
  /*
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
     */
  /*
     for (unsigned int i = 0; i < xics.size(); ++i) {
     ofs << xics[i].mz_ << "\t" << xics[i].rt_ << std::endl;
     }
     */
}
