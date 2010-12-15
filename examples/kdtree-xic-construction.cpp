/* $Id: example-xic-construction.cpp 1 2010-10-30 01:14:03Z mkirchner $
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

#include <array>
#include <iostream>
#include <fstream>
#include <tuple>
#include <sys/time.h>
#include <vector>

#include "spatial/kd_tree.h"

#include "fbi/connectedcomponents.h"

#include "example-xic-construction.h"
#include "example-xic-construction-opts.h"

int main(int argc, char* argv[])
{
  using namespace fbi;

  ProgramOptions options;
  if (!parseProgramOptions(argc, argv, options)) {
    return 0;
  }

  std::vector<Centroid> centroids = parseFile(options);

  timeval start, end;
  double mzWindowPpm = 2.0;
  double snWindow = 2.1;

  gettimeofday(&start, NULL);
  // construct kd-tree
  typedef std::array<double, 2> KeyType;
  typedef double MappedType;
  typedef ssrc::spatial::kd_tree<KeyType, MappedType> KdTree;
  KdTree kdtree;
  typedef std::vector<Centroid>::iterator CI;
  for (CI i = centroids.begin(); i != centroids.end(); ++i) {
    KeyType key;
    key[0] = i->mz_;
    key[1] = i->sn_;
    MappedType value = i->abundance_;
    kdtree[key] = value;
  }
  kdtree.optimize();
  typedef KdTree::const_iterator KCI;
  std::vector<std::set<size_t> > adjList(centroids.size());
  for (KCI i = kdtree.begin(); i != kdtree.end(); ++i) {
    // construct the range query
    KeyType llh = {{ i->first[0] * (1 - mzWindowPpm * 1E-6),
      i->first[1] - snWindow - 0.3 }};
    KeyType urh = {{ i->first[0] * (1 + mzWindowPpm * 1E-6),
      i->first[1] + snWindow + 0.3}};
    // run range query
    KCI b = kdtree.begin(llh, urh);
    // store connectivity
    for (KCI j = b; j != kdtree.end(); ++j) {
        size_t k = std::distance(b, i);
        size_t l = std::distance(b, j);
        adjList[k].insert(l);
        adjList[l].insert(k);
    }
  }
  kdtree.clear();
  gettimeofday(&end, NULL);
  std::cout << centroids.size() << "\t" << static_cast<double>(end.tv_sec - start.tv_sec) +
    static_cast<double>(end.tv_usec - start.tv_usec)* 1E-6 << std::endl;

  std::vector<long unsigned int> labels;
  findConnectedComponents(adjList, labels); 

  std::ofstream ofs(options.outputfileName_.c_str());
  ofs.setf(std::ios::fixed, std::ios::floatfield);
  for (size_t i = 0;i < centroids.size();++i) {
    Centroid& c = centroids[i];
    ofs << c.rt_ << "\t" << c.mz_ << "\t" << c.sn_ << "\t" 
      << c.abundance_ << "\t" << labels[i] << "\n";
  }
  return 0;
}


