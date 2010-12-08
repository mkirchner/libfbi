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

#include <iostream>
#include <fstream>
#include <tuple>
#include <sys/time.h>
#include <vector>

#include "fbi/tuplegenerator.h"
#include "fbi/fbi.h"
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

  gettimeofday(&start, NULL);
  auto results = SetA<Centroid, 1, 2>::
      intersect(centroids, BoxGenerator(2, 2.1), BoxGenerator(2, 2.1));
  gettimeofday(&end, NULL);
  std::cout << centroids.size() << "\t" << static_cast<double>(end.tv_sec - start.tv_sec) +
    static_cast<double>(end.tv_usec - start.tv_usec)* 1E-6 << std::endl;

  typedef SetA<Centroid, 1, 2>::IntType LabelType;
  std::vector<LabelType> labels;
  findConnectedComponents(results, labels); 

  std::ofstream ofs(options.outputfileName_.c_str());
  ofs.setf(std::ios::fixed, std::ios::floatfield);
  for (size_t i = 0;i < centroids.size();++i) {
    Centroid& c = centroids[i];
    ofs << c.rt_ << "\t" << c.mz_ << "\t" << c.sn_ << "\t" 
      << c.abundance_ << "\t" << labels[i] << "\n";
  }
  return 0;
}


