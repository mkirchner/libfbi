/* $Id: example-isotope-patterns.h 1 2010-10-30 01:14:03Z mkirchner $
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

#ifndef __LIBFBI_EXAMPLES_EXAMPLEISOTOPEPATTERNS_H__
#define __LIBFBI_EXAMPLES_EXAMPLEISOTOPEPATTERNS_H__

#include <utility>
#include <tuple>
#include "fbi/tuplegenerator.h"
#include "example-xic-construction.h"

struct Xic
{
  double rt_;
  double mz_;
  double abundance_;
  Xic() : rt_(0.0), mz_(0.0), abundance_(0.0) {}
};

namespace fbi {

template<>
struct Traits<Xic> : mpl::TraitsGenerator<double, double> {};

} //end namespace fbi

struct XicBoxGenerator
{
  template <size_t N>
  typename std::tuple_element<N, 
    typename fbi::Traits<Xic>::key_type>::type 
  get(const Xic&) const;

  double mzOffset_;
  double mzWindowPpm_;
  double rtOffset_;
  double rtWindow_;

  XicBoxGenerator(double mzWindowPpm, double rtWindow)
    : mzOffset_(0.0), mzWindowPpm_(mzWindowPpm),
      rtOffset_(0.0), rtWindow_(rtWindow)
  {}

  XicBoxGenerator(double mzOffset, double mzWindowPpm, 
    double rtOffset, double rtWindow)
    : mzOffset_(mzOffset), mzWindowPpm_(mzWindowPpm),
      rtOffset_(rtOffset), rtWindow_(rtWindow)
  {}
};

template <>
std::pair<double, double>  
XicBoxGenerator::get<0>(const Xic & centroid) const 
{
  return std::make_pair(
    mzOffset_ + centroid.mz_* (1 - mzWindowPpm_ * 1E-6), 
    mzOffset_ + centroid.mz_* (1 + mzWindowPpm_ * 1E-6));
}

template <>
std::pair<double, double>  
XicBoxGenerator::get<1>(const Xic & centroid) const
{
  return std::make_pair(
    rtOffset_ + centroid.rt_ - rtWindow_, 
    rtOffset_ + centroid.rt_ + rtWindow_);
}

#endif

