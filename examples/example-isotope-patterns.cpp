/* $Id: example-isotope-patterns.cpp 1 2010-10-30 01:14:03Z mkirchner $
 *
 * Copyright (c) 2010 Buote Xu <buote.xu@gmail.com>
 * Copyright (c) 2010 Marc Kirchner <marc.kirchner@childrens.harvard.edu>
 *
 * This file is part of libfbi.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without  restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions: 
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR  OTHER DEALINGS IN
 * THE SOFTWARE.
 */
#include <iostream>
#include <fstream>
#include <vector>

#include "boost/date_time/posix_time/posix_time.hpp"

#include "fbi/tuplegenerator.h"
#include "fbi/fbi.h"
#include "fbi/connectedcomponents.h"

#include "example-isotope-patterns.h"
#include "example-xic-construction-opts.h"

int main(int argc, char* argv[])
{
  using namespace fbi;
  using namespace boost::posix_time;

  ProgramOptions options;
  if (!parseProgramOptions(argc, argv, options)) {
    return 0;
  }

  std::vector<Centroid> centroids = parseFile(options);


  ptime start = microsec_clock::universal_time();
  SetA<Centroid,1,2>::ResultType centroidResults = SetA<Centroid, 1, 2>::
      intersect(centroids, BoxGenerator(2, 2.1), BoxGenerator(2, 2.1));
  ptime end = microsec_clock::universal_time();
  time_duration td = end - start;

  std::cout << "elapsed time in seconds: " 
    << td.total_seconds()
    << std::endl;

  std::cout << "finding connected components... ";
  typedef SetA<Centroid, 1, 2>::IntType LabelType;
  std::vector<LabelType> labels;
  LabelType nComponents = findConnectedComponents(centroidResults, labels); 
  std::cout << nComponents << " components found." << std::endl;

  typedef std::vector<std::vector<LabelType> > VecVec;
  VecVec groups(nComponents);
  typedef std::vector<LabelType>::iterator IT;
  LabelType k = 0;
  for (IT i = labels.begin(); i != labels.end(); ++i, ++k) {
    groups[*i-1].push_back(k);
  }
  std::vector<Xic> xics;
  typedef VecVec::iterator VIT;
  for (VIT i = groups.begin(); i != groups.end(); ++i) {
    Xic xic;
    for (IT j = i->begin(); j != i->end(); ++j) {
      xic.mz_ += centroids[*j].mz_;
      xic.rt_ += centroids[*j].rt_;
      xic.abundance_ += centroids[*j].abundance_;
    }
    double n = static_cast<double>(i->size());
    xic.mz_ /= n;
    xic.rt_ /= n;
    xics.push_back(xic);
  }

  std::vector<XicBoxGenerator> boxGenerators;
  for (double i = 0.0; i < 4; i = i + 1) {
    boxGenerators.push_back(XicBoxGenerator(i/3., 2, 0.0, 12.0));
  }



  // search for isotope pattern candidates
  typedef SetA<Xic, 0, 1> XicSet;
  
  start = microsec_clock::universal_time();
  
  XicSet::ResultType xicResults = XicSet::
    intersect(xics, XicBoxGenerator(0.0,2,0.0,12.0), boxGenerators);

  end = microsec_clock::universal_time();
  td = end - start;

  std::cout << "elapsed time in seconds: "
    << td.total_seconds()<< std::endl;

  std::cout << "finding connected components...";
  typedef XicSet::IntType IsotopeLabelType;
  std::vector<IsotopeLabelType> isotopeLabels;
  IsotopeLabelType nIsotopePatterns = findConnectedComponents(xicResults, isotopeLabels); 
  std::cout << nIsotopePatterns << " components found." << std::endl;

  std::ofstream ofs(options.outputfileName_.c_str());
  ofs.setf(std::ios::fixed, std::ios::floatfield);
  for (std::size_t i = 0;i < xics.size();++i) {
    Xic& c = xics[i];
    ofs << c.mz_ <<"\t" << c.mz_ * (1 - 2 * 1E-6) <<"\t" <<c.mz_ * (1 + 2* 1E-6) <<c.rt_ <<  "\t" 
      << c.abundance_ << "\t" << isotopeLabels[i] << "\n";
  }
  return 0;
}


