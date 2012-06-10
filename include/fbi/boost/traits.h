/* $Id: traits.h 1 2010-10-30 01:14:03Z mkirchner $
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

#ifndef __LIBFBI_INCLUDE_FBI_TRAITS_H__
#define __LIBFBI_INCLUDE_FBI_TRAITS_H__

#include <fbi/tuplegenerator.h>

namespace fbi {

/**
 * \class Traits
 * \brief A template for defining Traits of a specific Box class, has to define
 * its types / limits.
 * 
 * The user has to specialize this class template to his needs, as the
 * fbi::Tree will use the typedefs in the
 * \ref Traits class to get the correct types.
 * @verbatim
 * namespace fbi{
 *
 * template<>
 * struct Traits<CentroidType> :: mpl::TraitsGenerator<double, int>{};
 * \\^- shorthand for the following:
 *
 * template<>
 * struct Traits<CentroidType> 
 * {
 *  typedef std::tuple<std::pair<double, std::less<double> >,
 *    std::pair<int, std::less<int> >  > dim_type;
 *  typedef std::tuple<std::pair<double, double>, 
 *    std::pair<int, int> > key_type;
 *  static key_type getLimits()
 *  {
 *    return std::make_tuple(
 *      std::make_pair(std::numeric_limits<double>::min(), 
 *        std::numeric_limits<double>::max()),
 *      std::make_pair(std::numeric_limits<int>::min(),
 *        std::numeric_limits<int>::max())
 *    )
 *  }
 * enum {
 *   defined = 1
 * }
 * };
 * 
 * // This only works if std::numeric_limits is available for the chosen types,
 * // otherwise the user has to define them himself.
 * } //end namespace fbi 
 * @endverbatim
 * \note  
 *  \see \ref std::numeric_limits
 *  \see \ref mpl::TraitsGenerator
 *
 */


template<typename T> 
struct Traits: public mpl::TraitsGenerator<boost::mpl::void_>{
  enum {
    defined = 0
  };
}; 


} //end namespace fbi

#endif
