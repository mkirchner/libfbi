/* $Id: tuplegenerator.h 1 2010-10-30 01:14:03Z mkirchner $
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

#ifndef __LIBFBI_INCLUDE_FBI_TUPLEGENERATOR_H__ 
#define __LIBFBI_INCLUDE_FBI_TUPLEGENERATOR_H__ 

#include <utility> // for std::pair
#include <limits>
#include <vector>
//c++0x
#include <tuple> //for tuple_element and std::get
#include <array>

namespace fbi {

namespace mpl {

template <bool b>
struct Bool2Type {
  enum {value = b};
};

template <typename T>
struct Type2Type {
  typedef T OriginalType;
};


struct ConstructPair {
  template <typename T>
  struct apply { typedef std::pair<T,T>  type; };
};

struct ExtractFirst {
  template <typename T>
  struct apply { typedef typename T::first_type type; };
};

template <typename Metafun, typename Tuple>
struct mod;

template<typename Metafun, template<typename ...> class Tuple, 
  typename ... Types>
struct mod<Metafun, Tuple<Types...> > 
{
  typedef Tuple<typename Metafun::template apply<Types>::type...> type;
};

/** 
 * \brief Meta-template struct to add syntactic sugar for creating tuples of
 * pairs. 
 *
 * This is a meta-template way to convert a list of parameters into a tuple of
 * pairs, while also checking their type against the given AccessorType.  Let's
 * call the finished product EndTuple.
 * Idea: 
 * 1: Make typed tuple out of the parameter pack (the list). Let's call it TypeTuple.
 * It's only a typedef/template parameter and will allow us getting the type of
 * an individual component by using indices.
 * 2: By recursing, create a list of indices which correspond to the elements
 * which should be the heads of their respective pairs. (creating
 * (|parameter-pack| / 2) pairs.
 * 3: Now that we have the indices as a parameter pack, we can do the following:
 * 4: When expanding a parameter pack, it is possible to apply a function on
 * every single parameter by itself.
 * 5: Create a tuple from the parameter list.
 * 6: Now create a tuple by calling a "get" function on all of the
 * aforementioned tuple, one time for every index in the indices list, 
 * this get function should return a pair of items.
 * 6: Typechecking is done by trying to cast the resulting tuple into an
 * AccessorType, if they don't match the types were wrong.
 */


//That's a list of indices, as std::tuple can't handle std::size_t
template<std::size_t...>
struct Indices{};


/** 
 *\class PackImpl
 *Base class template, this class is calculating the indices list by recursing
 * over the TypeTuple while counting upwards (till I equals len(Tuple)/2
 * When it has finished recursing, the tail::index_values is equal to 
 * indices<0,2,4,...,len(Tuple) - 2>.
 * By defining the index_values as Next::index_values, the top index_values will
 * become equal to the tail index_values.
 *
 * The return type will be equal to a tuple of pairs, the pairs will have the
 * same first_type and second_type as the 0, 2, 4th ...element of the TypeTuple

*/
template<std::size_t I, class Indices, class Tuple, std::size_t N>
struct PackImpl;

template <std::size_t I, std::size_t... IndicesPack, class Tuple,
std::size_t N>
struct PackImpl<I, Indices<IndicesPack...>, Tuple, N>
{
  typedef PackImpl<I + 1, Indices<IndicesPack..., 2 * I>, Tuple, N>
      Next;
  typedef typename Next::index_values index_values;
  typedef typename Next::return_type return_type;

};

template <std::size_t N, std::size_t... IndicesPack, class Tuple>
struct PackImpl<N, Indices<IndicesPack...>, Tuple, N>
{
  typedef Indices<IndicesPack...> index_values;
  typedef std::tuple<
      std::pair<typename std::tuple_element<IndicesPack, Tuple>::type,
      typename std::tuple_element<IndicesPack, Tuple>::type
          >...
          > return_type;
};






template<std::size_t Index, class Tuple>
inline std::pair<typename std::tuple_element<Index, Tuple>::type,
       typename std::tuple_element<Index, Tuple>::type
       >
       make_pair(const Tuple& ts) {
         return 
             std::pair<
             typename std::tuple_element<Index, Tuple>::type,
         typename std::tuple_element<Index, Tuple>::type
             > (
                 std::get<Index>(ts), std::get<Index + 1>(ts));
       }

template<class Result, std::size_t... IndicesPack, class Tuple>
inline Result make_pair_tuple(Indices<IndicesPack...>, const Tuple& ts)
{
  return Result(make_pair<IndicesPack>(ts)...);

}

template<class... Types>
struct Pack {
  static_assert(((sizeof...(Types) & 1) == 0), "variadic parameter number must be even");
  typedef PackImpl<0, Indices<>, std::tuple<Types...>, sizeof...
      (Types)/2> Impl;
  typedef typename Impl::return_type return_type;
  typedef typename Impl::index_values index_values;
  static return_type Make(const Types&... ts) {
    return make_pair_tuple<return_type>(index_values(),
                                        std::make_tuple(ts...));
  }

};

/** 
 * \class TraitsGenerator
 * \brief A small generator struct to make writing the traits class easier.
 *  
 * \see fbi::Traits
 * @verbatim
 * namespace fbi{
 *
 * Traits<CentroidType> :: mpl::TraitsGenerator<double, double, int>{};
 *
 * };
 * // This only works if TraitsCreator is instantiated with numeric types, otherwise the user has to define them himself. //
 * } //end namespace fbi 
 * @endverbatim
 */
template <typename ... T>
struct TraitsGenerator{
  /**
   * Define the dimensions and their corresponding comparison, 
   * i.e. double, less<double>, int, less<int> for two dimensions
   */
  typedef std::tuple<std::pair<T, std::less<T> > ...> dim_type;
  /** 
   * As every dimension holds a pair of values, we define the correct tuples as key_type, 
   * using the types given in dim_type
   */
  typedef std::tuple<std::pair<T, T> ...> key_type;
  /**For every dimensionof the key_type, upper and lower bounds have to be defined */
  static key_type getLimits(){
    return std::make_tuple(std::make_pair(std::numeric_limits<T>::min(), std::numeric_limits<T>::max())...);
  }
    enum {
      defined = 1
    };
};


/** Base template class*/
template <size_t ...>
struct IndexChecker;

/** Check if any of the indices are bigger than T*/
template <size_t T, size_t FirstIndex, size_t ... Indices>
struct IndexChecker<T, FirstIndex, Indices...>{ 
  /** enum for compile-time checks*/
  enum {
    /** value == false if any index is greater or equal than T*/
    value = (FirstIndex < T) && IndexChecker<T, Indices...>::value
  };
};
/** Always return true if there are no indices left*/
template <size_t T>
struct IndexChecker<T>{
  /** enum for compile-time checks*/
  enum {
/** Always return true if there are no indices left*/
    value = true
  };
};


template <class TraitsType, std::size_t ... TIndices>
struct TypeExtractor {
  
  enum {
  /** TINDICESCORRECT is 1 when all TIndices can be used to access a dimension the 
    * key_type, otherwise it's 0
    */ 
    TINDICESCORRECT =
    IndexChecker<std::tuple_size<typename TraitsType::key_type>::
      value, TIndices...>::value && TraitsType::defined 
  };

  /** Empty TIndices shouldn't work */
  static_assert(
    sizeof...(TIndices) > 0, 
    "Please specify at least one index for the dimensions");

  /** Throw a compile error when the indices aren't correct. */
  static_assert(TINDICESCORRECT, 
                "Please check your SetA-Indices again, the \
                Traits<BoxType>::key_type \
                does not have enough dimensions to use your indices");

  template <bool Correct, class Placeholder>
  struct ExtractorImpl {
  
   typedef typename std::tuple < 
      typename std::tuple_element<
      TIndices, typename TraitsType::key_type 
      >::type ...
    > key_type;
  
  typedef typename std::tuple <
      typename std::tuple_element<
        TIndices, typename TraitsType::dim_type 
        >
      ::type::second_type ...
      > comp_type;   
      
  enum {
    ExtractionSuccessful = true
  };
  };

  template<class Placeholder>
  struct ExtractorImpl<false, Placeholder> {
    
    typedef std::array<void *, sizeof...(TIndices)> key_type;
    typedef std::array<std::pair<void, std::less<void *> >, sizeof...(TIndices) > comp_type;

  enum {
    ExtractionSuccessful = false
  };

  };
  typedef typename ExtractorImpl<TINDICESCORRECT, Type2Type<TraitsType> >::key_type key_type;
  typedef typename ExtractorImpl<TINDICESCORRECT, Type2Type<TraitsType> >::comp_type comp_type;

  enum {
    ExtractionSuccessful = ExtractorImpl<TINDICESCORRECT, Type2Type<TraitsType> >::ExtractionSuccessful
  };

};


/** 
 * Templated helper class, returns the number of functors passed.
 * As we should be able to pass multiple types of functors and even
 * multiple objects per type, we need to calculate how many functors there are.
 * This is needed to calculate the original indices of the boxes from given pointers,
 * as only they are passed around - the re-identification is made by pointer-arithmetic.
 */
struct FunctorChecker {

  /** 
   *If the current head is a single functor-object, add 1
   * \param functor A single functor, having a
   * \verbatim get<Dim>(const BoxType & ) const; \endverbatim method
   * \param functors Other functors, see functor
   */
  template <class Functor, class ...Functors>
  static std::size_t count(const Functor & functor, const Functors& ...functors) {
    return 1 + count(functors...);
  }

  /** If there is only a single functor left, add 1
   * \param functor A single functor, having a
   * \verbatim get<Dim>(const BoxType & ) const; \endverbatim method
   */
  template <class Functor>
  static std::size_t count(const Functor & functor) {
    return 1;
  }

  /** Specialized function call, if the passed object is a std::vector of functors,
   * count all of them as we will create one object per vector-size.
   * \param functor A single functor or a vector of functors of the same type,
   *  having a \verbatim get<Dim>(const BoxType & ) const; \endverbatim method
   * \param functors More functors
   * 
   */
  template <class Functor, class ...Functors>
  static
  std::size_t count(const std::vector<Functor> & functor, const Functors& ...functors) {
    return functor.size() + count(functors...);
  }
  /** 
   * Specialized case, if only one functor-type is left and that one is a vector,
   * count every instantiated object of that typein the container.
   * \param functor A single functor or a vector of functors of the same type,
   *  having a \verbatim get<Dim>(const BoxType & ) const; \endverbatim method
   *
   */
  template <class Functor>
  static std::size_t count(const std::vector<Functor> & functor) {
    return functor.size();
  }




};










} //end namespace mpl
} //end namespace libfbi


#endif //include guard
