/* $Id: simple-example.cpp 2 2010-10-30 01:16:54Z mkirchner $
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

#include <string>
#include "fbi/tuplegenerator.h" //TraitsGenerator
#include "fbi/fbi.h" //SetA::intersect

struct Location {
    double xPos;
    double yPos;
    std::string name;
    unsigned int nVisits;
};

std::vector<Location> locations;

namespace fbi {
  
template<>
struct Traits<Location> : mpl::TraitsGenerator<double, double> {};

} //end namespace

struct LocationBoxGenerator
{
  template <size_t N>
  typename std::tuple_element<N, 
    typename fbi::Traits<Location>::key_type>::type 
  get(const Location&) const;

  double xWidth_;
  double yWidth_;

  LocationBoxGenerator(double x, double y)
    : xWidth_(x), yWidth_(y) {}
};

template <>
std::pair<double, double>  
LocationBoxGenerator::get<0>(const Location& loc) const 
{
  return std::make_pair(
    loc.xPos - xWidth_ / 2, loc.xPos + xWidth_ / 2);
}

template <>
std::pair<double, double>  
LocationBoxGenerator::get<1>(const Location& loc) const 
{
  return std::make_pair(
    loc.yPos - yWidth_ / 2, loc.yPos + yWidth_ / 2);
}

int main() {
  auto adjList = fbi::SetA<Location, 0, 1>::intersect(
    locations, LocationBoxGenerator(0.5, 0.5), LocationBoxGenerator(0.5, 0.5));
}
