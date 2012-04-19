/* $Id: example-xic-construction.h 1 2010-10-30 01:14:03Z mkirchner $
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

#ifndef __LIBFBI_EXAMPLES_EXAMPLEXICCONSTRUCTION_H__
#define __LIBFBI_EXAMPLES_EXAMPLEXICCONSTRUCTION_H__

#include <utility>
#include <tuple>
#include "fbi/fbi.h"
#include "fbi/tuplegenerator.h"

struct ProgramOptions 
{
  double mzWindowLow_;
  double snWindowLow_;
  double mzWindowHigh_;
  double snWindowHigh_;
  std::string inputfileName_;
  std::string outputfileName_;
};


struct Centroid 
{
  double rt_;
  double mz_;
  double sn_;
  double abundance_;
  Centroid(const double &rt, const double& mz, const double & sn, 
    const double & abundance) 
    : rt_(rt), mz_(mz), sn_(sn), abundance_(abundance) {}
};

namespace fbi {

template<>
struct Traits<Centroid> : mpl::TraitsGenerator<
  double, double, double, double> {};

} //end namespace fbi

struct BoxGenerator
{
  template <size_t N>
  typename boost::tuples::element<N, 
    typename fbi::Traits<Centroid>::key_type>::type 
  get(const Centroid &) const;

  double mzOffset_;
  double mzWindowPpm_;
  double snOffset_;
  double snWindow_;

  BoxGenerator(double mzWindowPpm, double snWindow)
    : mzOffset_(0.0), mzWindowPpm_(mzWindowPpm),
      snOffset_(0.0), snWindow_(snWindow)
  {}

  BoxGenerator(double mzOffset, double mzWindowPpm, 
    double snOffset, double snWindow)
    : mzOffset_(mzOffset), mzWindowPpm_(mzWindowPpm),
      snOffset_(snOffset), snWindow_(snWindow)
  {}
};

template <>
std::pair<double, double>  
BoxGenerator::get<1>(const Centroid & centroid) const 
{
  return std::make_pair(
    mzOffset_ + centroid.mz_* (1-  mzWindowPpm_* 1E-6), 
    mzOffset_ + centroid.mz_* (1+ mzWindowPpm_ * 1E-6) );
}

template <>
std::pair<double, double>  
BoxGenerator::get<2>(const Centroid & centroid) const
{
  return std::make_pair(
    snOffset_ + centroid.sn_ - snWindow_ - 0.3, 
    snOffset_ + centroid.sn_ + snWindow_ + 0.3 );
}

std::vector<Centroid> parseFile(ProgramOptions & options)
{
  std::vector<Centroid> centroids;
  std::ifstream ifs(options.inputfileName_.c_str());
  ifs.setf(std::ios::fixed, std::ios::floatfield);

  double mz, rt, sn, abundance;

  while (ifs >> rt >> mz >> sn >> abundance)
  {
    if (mz >=options.mzWindowLow_ && mz <= options.mzWindowHigh_ &&
      sn >= options.snWindowLow_ && sn <= options.snWindowHigh_) {
        centroids.insert(centroids.end(), Centroid(rt, mz, sn, abundance));
    }
  }
  return centroids;
}








#endif

