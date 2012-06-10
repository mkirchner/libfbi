/* $Id$
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
 #ifndef __LIBFBI_INCLUDE_FBI_TUPLE_H__
#define __LIBFBI_INCLUDE_FBI_TUPLE_H__

#include <fbi/config.h>
#if defined(__USE_VARIADIC_TEMPLATES__)
#include <tuple>
namespace fbi{
template <std::size_t __i, typename _Tp, class Safety = std::pair<void*, void *> >
struct tuple_element {
 
enum {
VALID = (std::tuple_size<_Tp>::value > __i)
};
template <bool, typename>
struct TupleImpl {
  typedef typename std::tuple_element<__i, _Tp>::type type;
};

template <typename DUMMY>
  struct TupleImpl<false, DUMMY> {
  typedef Safety type;
  static_assert(sizeof(DUMMY) == 0, 
  "fbi::tuple_element is trying to access an index past the size of the tuple");
};

typedef typename TupleImpl<VALID, _Tp>::type type;

};

} //end namespace fbi
#else 
#include <boost/tuple/tuple.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/int.hpp>

namespace fbi{

template <bool, int N, class T, class Safety>
struct TupleImpl {
  typedef typename boost::tuples::element<N,T>::type type;
};

template <int N, class T, class Safety>
  struct TupleImpl<false, N, T, Safety> {
  enum {
    CHECK = sizeof(Safety) == 0
  };

  BOOST_MPL_ASSERT_MSG(CHECK, ACCESSING_INDEX_PAST_TUPLE, (boost::mpl::int_<N>, T));
  
  typedef Safety type;
};

template <int N, class T, class Safety = std::pair<void*, void *> >
struct tuple_element {
  enum {
    VALID = (boost::tuples::length<T>::value > N)
  };

  typedef typename TupleImpl<VALID, N, T, Safety>::type type;

};
} //end namespace fbi

#endif 


#endif
