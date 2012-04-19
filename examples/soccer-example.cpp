/* $Id: soccer-example.cpp 2 2010-10-30 01:16:54Z mkirchner $
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
#include <vector>
#include "fbi/tuple.h"
#include "fbi/tuplegenerator.h" //TraitsGenerator
#include "fbi/fbi.h" //SetA::intersect

struct Player {
    std::string name;
    unsigned int position;
    std::pair<double, double> price;
    unsigned int age;
};

struct Club {
    std::string name;
    double minBudget;
    double maxBudget;
};

namespace fbi {
  
template<>
struct Traits<Player> : mpl::TraitsGenerator<double, double, unsigned int> {};

template<>
struct Traits<Club> : mpl::TraitsGenerator<double, unsigned int, double> {};

} // end namespace

struct PlayerBoxGenerator
{
  template <std::size_t N>
  typename fbi::tuple_element<N, 
    typename fbi::Traits<Player>::key_type>::type 
  get(const Player&) const;
};

template <>
std::pair<double, double>  
PlayerBoxGenerator::get<0>(const Player& p) const 
{
  return std::make_pair(
    static_cast<double>(p.position), 
    static_cast<double>(p.position));
}

template <>
std::pair<double, double>  
PlayerBoxGenerator::get<1>(const Player& p) const 
{
  return p.price;
}

template <>
std::pair<unsigned int, unsigned int>  
PlayerBoxGenerator::get<2>(const Player& p) const 
{
  return std::make_pair(p.age, p.age);
}

class ClubBoxGenerator
{
  public:
  ClubBoxGenerator(const unsigned int position,
    const std::pair<unsigned int, unsigned int>& age)
  : position_(position), age_(age)
  { 
    age_.second++;
  }

  template <std::size_t N>
  typename fbi::tuple_element<N,
    typename fbi::Traits<Club>::key_type>::type
  get(const Club&) const;

  private:
    unsigned int position_;
    std::pair<unsigned int, unsigned int> age_;
};

template <>
std::pair<double, double>
ClubBoxGenerator::get<0>(const Club& c) const
{
    return std::make_pair(
      static_cast<double>(position_) - 0.1,
      static_cast<double>(position_) + 0.1);
}

template <>
std::pair<unsigned int, unsigned int>  
ClubBoxGenerator::get<1>(const Club& c) const 
{
  return age_;
}

template <>
std::pair<double, double>
ClubBoxGenerator::get<2>(const Club& c) const
{
    return std::make_pair(c.minBudget, c.maxBudget);
}

int main() {
  std::vector<Player> players;
  std::vector<Club> clubs;
  std::vector<ClubBoxGenerator> cbg;
  cbg.push_back(ClubBoxGenerator(7, std::make_pair(19, 24)));
  cbg.push_back(ClubBoxGenerator(11, std::make_pair(19, 24)));
  fbi::SetA<Player, 0, 1, 2>::ResultType adjList = fbi::SetA<Player, 0, 1, 2>::SetB<Club, 0, 2, 1>::intersect(
    players, PlayerBoxGenerator(), clubs, cbg);
  return 0;
}
