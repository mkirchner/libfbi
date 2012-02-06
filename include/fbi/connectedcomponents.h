/* $Id: connectedcomponents.h 1 2010-10-30 01:14:03Z mkirchner $
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

#ifndef __LIBFBI_INCLUDE_FBI_CONNECTEDCOMPONENTS_H__
#define __LIBFBI_INCLUDE_FBI_CONNECTEDCOMPONENTS_H__

#include <vector>
#include <boost/preprocessor.hpp>

#if _MSC_VER && !__INTEL_COMPILER
  #define MSWORKAROUND 1
#else
  #define MSWORKAROUND 0
#endif

template < class Container >
BOOST_PP_IIF(MSWORKAROUND, typename Container::value_type::template value_type, typename Container::value_type::value_type)
findConnectedComponents(const Container & adjacencyList,
  std::vector<
BOOST_PP_IIF(MSWORKAROUND, typename Container::value_type::template value_type, typename Container::value_type::value_type)
  >& labels)
{
  typedef typename Container::value_type::const_iterator IT; 
  typedef typename Container::value_type::value_type T;
  labels.resize(adjacencyList.size(), 0);
  T currentLabel = 1;
  std::vector<T> stack;

  for (T nodeCounter = 0; 
    nodeCounter != static_cast<T>(adjacencyList.size()); ++nodeCounter) {
    if (labels[nodeCounter] != 0) 
      continue;
    stack.push_back(nodeCounter);
    while (!stack.empty()) {
      const T currentNode = stack.back();
      labels[currentNode] = currentLabel;
      stack.pop_back();
      IT it = adjacencyList[currentNode].begin();
      const IT end = adjacencyList[currentNode].end();
      for (; it != end; ++it) {
        if (labels[*it] == 0) stack.push_back(*it);
      }
    }
    ++currentLabel;
  }
  return currentLabel-1;
}
#undef MSWORKAROUND

#endif
