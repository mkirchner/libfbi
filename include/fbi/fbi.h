/* $Id: fbi.h 1 2010-10-30 01:14:03Z mkirchner $
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

#ifndef __LIBFBI_INCLUDE_FBI_HYBRID_H__
#define __LIBFBI_INCLUDE_FBI_HYBRID_H__

//C99
#include <stdint.h>
#include <stdlib.h>
//C++
#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <utility>
#include <vector>
#include <iostream>
//c++0x
#include <random>
#include <tuple>
//tree
#include <fbi/config.h>
#include <fbi/traits.h>

#include <fbi/tuplegenerator.h>

namespace fbi {

  /**
   * \class SetA
   *
   * \brief A class to find intersections between cartesian products of 
   *   arbitrarily-typed intervals.
   * 
   * Based on "Fast Software for Box Intersections", 
   *  by Afra Zomorodian, Herbert Edelsbrunner,
   * 
   * \tparam BoxType The objects we're looking at, 
   *  Traits<BoxType> has to available. 
   * \note To work correctly on the given types, 
   *  the templated Traits class has to be specialized
   *  for the given BoxType (and QBoxType, if needed).
   *
   * \see Traits 
   * \tparam TIndices Two boxes shall intersect if they do in these 
   *  dimensions, given as a parameter pack of std::size_t.
   * \note The user has to specify at least one index to work on, 
   *  empty TIndices won't return any results.
   */



template <typename BoxType, std::size_t ... TIndices>
class SetA{



 private:
  /** Keep ctor private to leave it as a class - 
    *   we don't want a tree object.
    */
  SetA(); 

  typedef SetA<BoxType, TIndices...> SETA;

  /** A compile-time constant to mark recursion tails.*/
  enum {LASTDIM = sizeof...(TIndices) - 1}; 

  /** Short typedef to offer the STL value_type identifier for this tree.*/
  typedef BoxType value_type;

  /** TINDICES is 1 when all TIndices can be used to access a dimension the 
    * key_type, otherwise it's 0
    */ 
  enum {
    TINDICESCORRECT =
    mpl::IndexChecker<std::tuple_size<typename Traits<BoxType>::key_type>::
      value, TIndices...>::value 
  };

  /** \brief Empty TIndices shouldn't work */
  static_assert(
    sizeof...(TIndices) > 0, 
    "Please specify at least one index");

  /** \brief Throw a compile error when the indices aren't correct. */
  static_assert(TINDICESCORRECT, 
                "Please check your SetAIndices again, the \
                Traits<BoxType>::key_type \
                does not have enough dimensions to use your indices");


  /**
   *  The key is represented by a pair of values for each dimension that 
   *  should be considered for intersection tests, 
   * these pairs are wrapped in a std::tuple for easy access via std::get.
   * \see \ref fbi::Traits
   * \note This key_type is a subset of the original key_type, 
   * depending on the indices we're using.
   */
  typedef typename std::tuple < 
      typename std::tuple_element<
      TIndices, typename Traits<value_type>::key_type 
      >::type ...
    > key_type;


  /**
   * As every dimension can have a different type, 
   * the appropriate comparison-operator
   * has to be known during compile-time. 
   * \see \ref fbi::Traits
   */
  typedef typename std::tuple <
      typename std::tuple_element<
        TIndices, typename Traits<value_type>::dim_type 
        >
      ::type::second_type ...
      > comp_type;      

 public:
  /** 
   * Every intersection between two elements will be represented by their 
   * corresponding index in the container of type C, 
   * to save memory we decided to artificially 
   * limit the amount of elements which allows us
   * to use 32Bit integers instead of 64Bit.
   *
   */
  typedef uint32_t IntType; 

  /**
   * The intersections will be returned as a adjacency list, 
   * to prevent parallel edges we're using a set for the inner dimension.
   */

  typedef std::vector<std::set<IntType> > ResultType; 


  /** 
    * \class SetB
   * \brief Subclass, when the type of the query objects differs from the 
   * type of the data objects, we need an 
   * additional layer of template specialization.
   */

  template <typename QueryType, std::size_t ... QIndices>
  struct SetB;

 private:

  /** 
    * \class KeyCreator
   * \brief Handle the creation of keys by using given functors. 
   * As we have to handle indices which aren't 
   * known at the time of the tree instantiation, we have to add 
   * another template class to use different std::size_t parameter packs.
   *
   */
  template <std::size_t ...KeyCreatorIndices>
  struct KeyCreator;


  /** 
  * \class State
   * \brief Holding various information like the random generator state, 
   * some specific addresses of data and dimensional limits.
   * 
   */
  class State;

  /**
  * \class HybridScanner
   * \brief Find intersections between boxes by 
   * recursing through virtual segment trees.
   *
   * We're interested in finding the intersections between 
   * two bipartite sets.
   * Similar to other divide & conquer algorithms, 
   * we will recursively divide the problem into smaller pieces, 
   * till we only have to check for intersections between 
   * two small subsets with a brute-force algorithm.
   * To find all intersections, we have to alternate between 
   * looking at the intervals as points and as intervals,
   * which means that every instance of HybridScanner has to spawn 
   * two HybridScanners for the next dimension with
   * alternating inputs.
   * 
   * \note Like Quicksort, we would like to use the median to 
   * divide our sets into equal parts.
   * As an exact median calculation would be too tiresome, we used a 
   * heuristic which is a generalized version of the
   * standard medianOfThree method to find a good approximation. 
   * \see ref getApproxMedian
   * \see ref OneWayScanner
   */
  template <bool PointsContainQueries, std::size_t Dim>
  struct HybridScanner;


  /**
    * \class OneWayScanner
   * \brief A scanner to find one type of intersection between 
   * two sets of intervals by looking at one of them as a set of points. 
   * 
   * Looking at the one-dimensional problem, two intervals can 
   * only overlap if at least one lower endpoint
   * (also called the Head) is contained inside the other interval.
   * We're finding intersections by iterating through the 
   * sorted set of points while looking at possible intervals
   * which could contain a given point.
   * The complexity of the problem has, given two sets of size n and m, 
   * an upper bound of \f$ O(n * m) \f$, 
   * for the general case we are looking at
   * \f$ O(m * log(m/n) + m + n ) \f$.
   * 
   * \note As we're only looking at one combination of 
   * points/intervals here during one call, only one
   * half of the possible combinations is found, 
   * this is why OneWayScanner should always be in pairs.
   * \tparam PointsContainQueries As we have to alternate between views, 
   * we would like to know which of the two sets at the beginning
   * our points are a subset of. 
   * This is needed to calculate the resulting indices by 
   * subtracting pointers.
   * \tparam Dimension where we should first look for an intersection, 
   * all following dimensions will be 
   * evaluated via comparisons between two objects.
   */
  template <bool PointsContainQueries, std::size_t Dim>
  struct OneWayScanner;

  /** 
   * \class lessHead
   * \brief This Functor compares the lower endpoints of 
   * two key_type objects in a given dimension.
   * \see std::pair<T>::operator<() 
   * \tparam Dim dimension to compare in.
   */ 
  template <std::size_t Dim>
  struct lessHead;

  /** 
   * \class lessTail
   * \brief This Functor compares the upper endpoints of 
   * two key_type objects in a given dimension.
   * \see std::pair<T>::operator<() 
   * \tparam Dim dimension to compare in.
   */

  template <std::size_t Dim>
  struct lessTail;

  /** 
   * \class IntersectionTester
   * \brief Check if two objects intersect in 
   * all dimensions between (including) Dim and Limit.
   *
   * If the amount of data is small enough, 
   * straightforward checking for two objects can be used.
   * Note that this should only be used after 
   * sufficient filtering of possible intersection candidates by the scanners.
   *
   */

  template <std::size_t Dim, std::size_t Limit>
  struct IntersectionTester;


  /** 
   * \class KeyPrinter 
   * \brief For Debug reasons, print all dimensions of a given key via cout.
   */
  template <std::size_t Dim, std::size_t Limit>
  struct KeyPrinter;

 public:



  /**
   * \brief Create two sets of keys by intersecting two sets of functors on a 
   * given dataContainer, and check if there are intersections 
   * between these keys in the dimensions the tree was initialized with.
   *
   * In the general case, we would like to have 
   * two different types of boxes, which, given functors, 
   * produce the same key_type. As their keys are then equal and also a 
   * representation of hyperrectangles 
   * (cartesian products of intervals), it is 
   * possible to test for intersections between them.
   * 
   * \param dataContainer A container with a const_iterator, 
   * it has to hold the same type as the one the tree was instantiated with.
   * Both key sets will be extracted by using 
   * functors on every value in this container.
   * \param ifunctor This has to be a class with a public 
   * \verbatim get<Dim>(const BoxType & ) const; \endverbatim method.
   * The tree will create the keys by extracting the indexed dimensions from the
   *  BoxType.
   * \param qfunctors Like ifunctor, 
   * these functors will each create a different key.
   * 
   * \return We'll return a std::vector<std::set<IntType> >, 
   * as parallel edges are possible but not desired for the end result.
   * An edge (pair of two boxes) consists of 2 IntType values, 
   * the head is implied by the index in the std::vector, as 
   * every vertex is holding a set of tail indices.
   *
   * \see \ref SetB::intersect
   * \callgraph
   */
  template <
  class BoxContainer,
        typename = 
          typename std::enable_if<
            std::is_same<typename BoxContainer::value_type, value_type>::value>
          ::type,
        typename IntervalFunctor,
        typename ... QueryFunctors
          >
          static
          ResultType intersect(
            const BoxContainer & dataContainer,
            const IntervalFunctor & ifunctor,
            const QueryFunctors & ... qfunctors
            )
          {
            return SetB<BoxType, TIndices...>::
                intersect(dataContainer, ifunctor, dataContainer, qfunctors...);
          }
/**
   * \brief Create two sets of keys by intersecting two sets of functors on a 
   * given dataContainer, and check if there are intersections 
   * between these keys in the dimensions the tree was initialized with.
   *
   * In the general case, we would like to have 
   * two different types of boxes, which, given functors, 
   * produce the same key_type. As their keys are then equal and also a 
   * representation of hyperrectangles 
   * (cartesian products of intervals), it is 
   * possible to test for intersections between them.
   * 
   * \param dataContainer A container with a const_iterator, 
   * it has to hold the same type as the one the tree was instantiated with.
   * Both key sets will be extracted by using 
   * functors on every value in this container.
   * \param ifunctor This has to be a class with a public 
   * \verbatim get<Dim>(const BoxType & ) const; \endverbatim method.
   * The tree will create the keys by extracting the indexed dimensions from the
   *  BoxType.
   * \param qfunctors Like ifunctor, 
   * these functors will each create a different key.
   * 
   * \return We'll return a std::vector<std::set<IntType> >, 
   * as parallel edges are possible but not desired for the end result.
   * An edge (pair of two boxes) consists of 2 IntType values, 
   * the head is implied by the index in the std::vector, as 
   * every vertex is holding a set of tail indices.
   *
   * \see \ref SetB::intersect
   * \callgraph
   */
  template <
  class BoxContainer,
        typename = 
          typename std::enable_if<
            std::is_same<typename BoxContainer::value_type, value_type>::value>
          ::type,
        typename IntervalFunctor,
        typename ... QueryFunctors
          >
          static
          ResultType thetaIntersect(
            const std::size_t cutoff,
            const BoxContainer & dataContainer,
            const IntervalFunctor & ifunctor,
            const QueryFunctors & ... qfunctors
            )
          {
            return SetB<BoxType, TIndices...>::
                thetaIntersect(cutoff, dataContainer, ifunctor, dataContainer, qfunctors...);
          }




 private:

  /**
   * Calculate an approximate median by using the median-of-three method
   * on a tertiary tree of given height.
   *
   * \tparam Dim Get the median of a vector of values in the given dimension.
   * \see \ref heuristicHeight()
   * \param[in] container The vector of values we're interested in
   * \param[in] height The depth of the triadic tree we're 
   * resolving to get a median
   * \param[in,out] state The object holding the random number generator
   * \param[in] less The comparison operator so that we can 
   * compare in a specific dimension
   *
   */
  template <std::size_t Dim>
  const 
  static 
  typename std::tuple_element<Dim, key_type>::type::first_type
  getApproxMedian(
    const std::vector<const key_type *> & container,
    const std::size_t height, State & state, 
    const typename std::tuple_element<Dim, comp_type>::type & less)
  {
    typedef typename std::tuple_element<Dim, key_type>::type::first_type 
      ValType;
    
    if (height == 0) {
      std::size_t random = state.randInt(0, container.size()-1);
      std::size_t rand1 = state.randInt(0,1);
      if (rand1 == 0)
      {
        return getHead<Dim>(container[random]);
      }
      else
      {
        return getTail<Dim>(container[random]);
      }
    }
    ValType t1 = getApproxMedian<Dim>(container, height-1, state, less);
    ValType t2 = getApproxMedian<Dim>(container, height-1, state, less);
    ValType t3 = getApproxMedian<Dim>(container, height-1, state, less);
    return medianOfThree(t1, t2, t3, less);
  }


  /** 
   * Create a vector with pointers to our interval data to save memory.
   *
   * \param container A STL container holding key_type objects, 
   * has to have a const_iterator.
   * As our \ref ScanFunctions have to pass and copy containers of intervals, 
   * using a pointer based approach was deemed 
   * faster and safer as we won't modify the actual data.
   *
   * \see \ref createIntervalVector()
   * \note Don't move/insert into the intervalContainer afterwards, as it will 
   * invalidate all pointers.
   */
  
  static const std::vector<const key_type *>
  createPtrVector (const std::vector<key_type> & container) {
    std::vector<const key_type *> ptrVector(container.size());
    typename std::vector<key_type>::const_iterator it = container.begin();
    std::size_t i = 0;
    while (it != container.end()){
      ptrVector[i] = &(*it);
      ++it;
      ++i;
    }
    return ptrVector;
  }


  /**
   *
   * Some easy ways to get the correct interval / endpoints in a given dimension
   * \tparam Dim dimension to work on.
   */

  /** Comfort function to get the interval in the correct dimension */
  template <std::size_t Dim>
  static inline typename std::tuple_element<Dim,key_type>::type 
  getKey(const key_type & key) {
    return std::get<Dim>(key);
  }

  /** Comfort function to get the interval in the correct dimension */
  template <std::size_t Dim>
  static inline typename std::tuple_element<Dim,key_type>::type 
  getKey(const key_type * key) {
    return std::get<Dim>(*key);
  }

  /** Comfort function to get the lower endpoint in the correct dimension */
  template <std::size_t Dim>
  static inline typename std::tuple_element<Dim,key_type>::type::first_type 
  getHead(const key_type & key) {
    return getKey<Dim>(key).first;
  }

  /** Comfort function to get the lower endpoint in the correct dimension */
  template <std::size_t Dim>
  static inline typename std::tuple_element<Dim,key_type>::type::first_type
  getHead(const key_type * key) {
    return getKey<Dim>(key).first;
  }

  /** Comfort function to get the upper endpoint in the correct dimension */
  template <std::size_t Dim>
  static inline typename std::tuple_element<Dim,key_type>::type::first_type
  getTail(const key_type & key) {
    return getKey<Dim>(key).second;
  }



  /** Comfort function to get the upper endpoint in the correct dimension */
  template <std::size_t Dim>
  static inline typename std::tuple_element<Dim,key_type>::type::first_type
  getTail(const key_type * key){
    return getKey<Dim>(key).second;
  }


  /** Comfort function to get the comparison functor in the correct dimension */
  template <std::size_t Dim>
  static inline const typename std::tuple_element<Dim, comp_type >::type
  getCompareFunctor(){
    return typename std::tuple_element<Dim,comp_type>::type();
  }


  /** Sort a container of pointers to key_type objects, compare their 
    * lower endpoints in the specified dimension 
    */ 
  template <std::size_t Dim>
  static inline void
  sortContainerHead(std::vector<const key_type * > & container){
    std::sort(container.begin(), container.end(), lessHead<Dim>());
  }

  /** Sort a container of pointers to key_type objects, compare their 
    * upper endpoints in the specified dimension 
    */ 
  template <std::size_t Dim>
  static inline void
  sortContainerTail( std::vector<const key_type * > & container){
    std::sort(container.begin(), container.end(), lessTail<Dim>());
  }
  /**
   * Calculate the median of three values, comparison functor has to be
   * provided.
   *
   * \tparam T type of compared values
   * \tparam Comp comparison functor
   *  \see \ref getApproxMedian() 
   *
   */
  template <typename T, typename Comp>
  static const T & medianOfThree (
      const T & t1, 
      const T & t2, 
      const T & t3, 
      const Comp & less ) 
  {
    if (less(t2, t1))
    {
      if (less(t3,t2)) return t2;
      if (less(t1,t3)) return t1;
      return t3;  
    }
    if (less(t3,t1)) return t1;
    if (less(t2,t3)) return t2;

    return t3; 

  }

  /**
   * Calculate the depth of the approximate median tree, as some 
   * constants here have to be found by experiment.
   *
   * \see \ref getApproxMedian() 
   */
  static const std::size_t inline heuristicHeight(std::size_t numElements)
  {
    return (size_t)(log((double)numElements));
  }

}; //end class tree





template <typename BoxType, std::size_t ... TIndices>
template <typename QBoxType, std::size_t ... QIndices>
struct SetA<BoxType, TIndices...>::
SetB {
 private:
  static_assert(sizeof...(QIndices) == sizeof...(TIndices), 
    "Your number of query-dimensions doesn't match your initial indices");

  typedef QBoxType qvalue_type;
  typedef typename std::tuple <
      typename std::tuple_element<
      QIndices, typename Traits<qvalue_type>::key_type
      >::type ...
      > qkey_type;
  static_assert(std::is_same<key_type, qkey_type>::value, 
    "Keytypes don't match");


  typedef typename std::tuple <
        typename std::tuple_element<QIndices, typename Traits<qvalue_type>::dim_type
      >::type::second_type ...
    > qcomp_type;      
  static_assert(std::is_same<comp_type, qcomp_type>::value, 
    "CompTypes don't match");

 public:
  /**
   * \callgraph
   * \brief The public interface for the user to start the algorithm. 
   *
   * This function initializes the keys (consisting of intervals in 
   * several dimensions) of the intervalSet and the querySet, along with 
   * pointers to pass them to separation algorithms.
   * The IndexCalculator has to be used to 
   * recalculate the original indices from the pointers.
   *
   *
   * \param[in] dataContainer STL Container providing a 
   *  forward iterator and holding a value_type.
   * \param[in] ifunctor A functor to operate on value_type and create a 
   *  key_type by getting a pair in every indexed dimension.
   * \param[in] qdataContainer STL Container providing a 
   *  forward iterator and holding a qvalue_type.
   * \param[in] qfunctors a variadic list of functors to create several 
   *  key_type objects which will be associated with the original data.
   * \return An adjacency list implemented by a random-accessable 
   *  outer container representing the vertices, 
   *  each holding a container of non-parallel edges.
   */

  template <
  class BoxContainer,
        typename = typename std::enable_if<std::is_same<typename BoxContainer::value_type, value_type>::value>::type,
        class QContainer,
        typename = typename std::enable_if<std::is_same<typename QContainer::value_type, qvalue_type>::value>::type,
        typename IntervalFunctor, 
        typename ... QueryFunctors
  > static
  ResultType intersect(
      const BoxContainer & dataContainer, 
      const IntervalFunctor & ifunctor, 
      const QContainer & qdataContainer,
      const QueryFunctors& ... qfunctors
      ) {
    static_assert( (sizeof...(QueryFunctors) > 0), 
      "Need at least one query functor.");
    if (dataContainer.empty()) { return ResultType();}
    // Generate the set of query boxes. The BoxType is an arbitrary,
    // user-specified type, that does not necessarily have any notion of
    // dimensionality. This call converts the BoxType data into the 
    // K-dimenstional boxes for fast box intersection.
    const auto queryIntervalVector = KeyCreator<QIndices...>::
      getVector(qdataContainer, qfunctors...);
    // Generate the set of data boxes. See above, just for the QueryBoxType.
    const auto dataIntervalVector = KeyCreator<TIndices...>::
      getVector(dataContainer, ifunctor);

    key_type limits = 
      make_tuple(
        std::get<TIndices>(Traits<value_type>::getLimits())
      ... );

    //const std::size_t numQueryFunctors = sizeof...(QueryFunctors); 
    const std::size_t numQueryFunctors = 
        mpl::FunctorChecker::count(qfunctors...); 
    //if we're looking at two different sets, use different indices for the elements!
    const std::size_t offset = 
        (reinterpret_cast<const char* const>(&(dataContainer)) == 
        reinterpret_cast<const char* const>(&(qdataContainer))) ? 0 : dataContainer.size();
    State state(
        limits,
        numQueryFunctors,
        &(queryIntervalVector[0]),
        &(dataIntervalVector[0]),
        offset
        );

    // Create a vector of pointers that reference the above query boxes. This
    // allows us to work on pointers and save a bit of memory.
    std::vector<const key_type *> pointsPtrVector = 
      createPtrVector(queryIntervalVector);
    std::vector<const key_type *> intervalsPtrVector = 
      createPtrVector(dataIntervalVector);

    ResultType resultVector(offset + qdataContainer.size());

    auto dimLimits = std::get<0>(state.getLimits()); 

    // Call the hybrid algorithm for stabbing queries in the interval vector.
    HybridScanner<true, 0>::
      scan(
        pointsPtrVector, 
        intervalsPtrVector, 
        dimLimits.first, 
        dimLimits.second,state, 
        resultVector 
      );
    // Reverse the previous call: queries in the "point" vector.
    HybridScanner<false,0>::
      scan(
        intervalsPtrVector, 
        pointsPtrVector, 
        dimLimits.first, 
        dimLimits.second,state, 
        resultVector
      );
    return resultVector;
  }
   /**
   * \callgraph
   * \brief The public interface for the user to start the algorithm. 
   *
   * This function initializes the keys (consisting of intervals in 
   * several dimensions) of the intervalSet and the querySet, along with 
   * pointers to pass them to separation algorithms.
   * The IndexCalculator has to be used to 
   * recalculate the original indices from the pointers.
   *
   * \param[in] cutoff The theta cutoff value for switching into OneWayScan
   * \param[in] dataContainer STL Container providing a 
   *  forward iterator and holding a value_type.
   * \param[in] ifunctor A functor to operate on value_type and create a 
   *  key_type by getting a pair in every indexed dimension.
   * \param[in] qdataContainer STL Container providing a 
   *  forward iterator and holding a qvalue_type.
   * \param[in] qfunctors a variadic list of functors to create several 
   *  key_type objects which will be associated with the original data.
   * \return An adjacency list implemented by a random-accessable 
   *  outer container representing the vertices, 
   *  each holding a container of non-parallel edges.
   */
  template <
  class BoxContainer,
        typename = typename std::enable_if<std::is_same<typename BoxContainer::value_type, value_type>::value>::type,
        class QContainer,
        typename = typename std::enable_if<std::is_same<typename QContainer::value_type, qvalue_type>::value>::type,
        typename IntervalFunctor, 
        typename ... QueryFunctors
  > static
  ResultType thetaIntersect(
      const size_t cutoff,
      const BoxContainer & dataContainer, 
      const IntervalFunctor & ifunctor, 
      const QContainer & qdataContainer,
      const QueryFunctors& ... qfunctors
      ) {
    static_assert( (sizeof...(QueryFunctors) > 0), 
      "Need at least one query functor.");
    if (dataContainer.empty()) { return ResultType();}
    // Generate the set of query boxes. The BoxType is an arbitrary,
    // user-specified type, that does not necessarily have any notion of
    // dimensionality. This call converts the BoxType data into the 
    // K-dimenstional boxes for fast box intersection.
    const auto queryIntervalVector = KeyCreator<QIndices...>::
      getVector(qdataContainer, qfunctors...);
    // Generate the set of data boxes. See above, just for the QueryBoxType.
    const auto dataIntervalVector = KeyCreator<TIndices...>::
      getVector(dataContainer, ifunctor);

    key_type limits = 
      make_tuple(
        std::get<TIndices>(Traits<value_type>::getLimits())
      ... );

    //const std::size_t numQueryFunctors = sizeof...(QueryFunctors); 
    const std::size_t numQueryFunctors = 
        mpl::FunctorChecker::count(qfunctors...); 
    //if we're looking at two different sets, use different indices for the elements!
    const std::size_t offset = 
        (reinterpret_cast<const char* const>(&(dataContainer)) == 
        reinterpret_cast<const char* const>(&(qdataContainer))) ? 0 : dataContainer.size();
    State state(
        limits,
        numQueryFunctors,
        &(queryIntervalVector[0]),
        &(dataIntervalVector[0]),
        offset,
        cutoff
        );

    // Create a vector of pointers that reference the above query boxes. This
    // allows us to work on pointers and save a bit of memory.
    std::vector<const key_type *> pointsPtrVector = 
      createPtrVector(queryIntervalVector);
    std::vector<const key_type *> intervalsPtrVector = 
      createPtrVector(dataIntervalVector);

    ResultType resultVector(offset + qdataContainer.size());

    auto dimLimits = std::get<0>(state.getLimits()); 

    // Call the hybrid algorithm for stabbing queries in the interval vector.
    HybridScanner<true, 0>::
      scan(
        pointsPtrVector, 
        intervalsPtrVector, 
        dimLimits.first, 
        dimLimits.second,state, 
        resultVector 
      );
    // Reverse the previous call: queries in the "point" vector.
    HybridScanner<false,0>::
      scan(
        intervalsPtrVector, 
        pointsPtrVector, 
        dimLimits.first, 
        dimLimits.second,state, 
        resultVector
      );
    return resultVector;
  }
}; //end class SetB

template <typename BoxType, std::size_t ... TIndices>
template <std::size_t ... KeyCreatorIndices>
struct SetA<BoxType, TIndices...>::
KeyCreator{

  /**
   * Calculate the intervals every value_type object is representing and save it
   * as a key_type object to work with.
   * \callgraph
   * 
   * \param functors Functors which provide a 
   * std::pair<...> get<std::size_t>(value_type) member function.
   * We create multi-dimensional intervals by using the functor on all 
   * indexed dimensions.
   * \param container A STL container with value_type objects, 
   * has to provide a forward iterator.
   * \see \ref createKeyType()
   */

  template <class Container, class ... Functors>
  static std::vector<key_type>
  getVector(const Container & container, const Functors& ...functors){
    static_assert(sizeof...(Functors) > 0, 
      "You need at least one functor to access your objects"); 
    typename Container::const_iterator it = container.begin();
    std::vector<key_type> intervalVector(
        container.size()* mpl::FunctorChecker::count(functors...)); 

    typename std::vector<key_type>::iterator intervalIt= intervalVector.begin();
    while (it != container.end())
    {
      createKeys(intervalIt, *it, functors...);
      ++it;
    }
    return intervalVector;
  }

  /**
   * This is the tail of a recursive call which uses different 
   * functors to produce several queries.
   *
   * \param intervalIt Iterator to a container which will 
   *  hold the different keys.
   * \param dataValue That's the box representing the class we're working on.
   * \param functor Functor which provides a 
   *  std::pair<...> get<std::size_t>(value_type) member function
   *  \see \ref createKeyType()
   */

  template <typename T, typename Functor>
  static void
  createKeys(typename std::vector<key_type>::iterator & intervalIt,
             const T & dataValue, 
             const Functor & functor) {
    *intervalIt = createKey(dataValue, functor);
    ++intervalIt; //we have to increment it so that the while loop is ok.
  }

  /**
   *  Use functors on a given value_type object to 
   * create a key for every functor.
   *
   * \param intervalIt Iterator to a container which will 
   *  hold the different keys.
   * \param dataValue That's the box representing the class we're working on.
   * \param functor Functor which provides a 
   *  std::pair<...> get<std::size_t>(value_type) member function
   * \param functors Other functors which will be used recursively.
   * \see \ref createKeyType()
   */
  template <typename T, typename Functor, typename ...Functors>
  static void 
  createKeys(
      typename std::vector<key_type>::iterator & intervalIt,
      const T & dataValue, 
      const Functor & functor, 
      const Functors & ...functors) {
    *intervalIt = createKey(dataValue, functor); 
    ++intervalIt;
    createKeys(intervalIt, dataValue, functors...);
  }

  /**
   * This is the tail of a recursive call which uses different 
   * functors to produce several queries.
   *
   * \param intervalIt Iterator to a container which will 
   *  hold the different keys.
   * \param dataValue That's the box representing the class we're working on.
   * \param functor Functor which provides a 
   *  std::pair<...> get<std::size_t>(value_type) member function
   *  \see \ref createKeyType()
   */

  template <typename T, typename Functor>
  static void
  createKeys(typename std::vector<key_type>::iterator & intervalIt,
             const T & dataValue, 
             const std::vector<Functor> & functor) 
  {
    for (std::size_t i = 0; i < functor.size(); ++i) {
      *intervalIt = createKey(dataValue, functor[i]);
      ++intervalIt; //we have to increment it so that the while loop is ok.
    }
  }

  /**
   *  Use functors on a given value_type object to 
   * create a key for every functor.
   *
   * \param intervalIt Iterator to a container which will 
   *  hold the different keys.
   * \param dataValue That's the box representing the class we're working on.
   * \param functor Functor which provides a 
   *  std::pair<...> get<std::size_t>(value_type) member function
   * \param functors Other functors which will be used recursively.
   * \see \ref createKeyType()
   */
  template <typename T, typename Functor, typename ...Functors>
  static void 
  createKeys(
      typename std::vector<key_type>::iterator & intervalIt,
      const T & dataValue, 
      const Functor & functor, 
      const std::vector<Functors> & ...functors) 
  {
    for (std::size_t i = 0; i < functor.size(); ++i) {
      *intervalIt = createKey(dataValue, functor[i]);
      ++intervalIt; //we have to increment it so that the while loop is ok.
    }
    createKeys(intervalIt, dataValue, functors...);
  }

  /**
   * \brief Create a key_type object (which is a tuple of pairs) by calling the 
   *  getInterval function of the box object in all relevant dimensions.
   *
   * \see \ref ScanFunctions
   */
  template <typename T, typename Functor>
  static inline const
  key_type 
  createKey(const T& val, const Functor & functor) {
    return std::make_tuple(functor.template get<KeyCreatorIndices>(val)...);
  }

}; //end struct KeyCreator




template <typename BoxType, std::size_t ... TIndices>
class SetA<BoxType, TIndices...>::
State
{
 private: 
  static std::size_t defaultHeightCalculator_(const std::size_t n) {
    return 
      static_cast<std::size_t> (
        std::max( 0., -3.0 + 1.8 * log10( static_cast<double>( n ) ) )
      ); 
  }

  const key_type limits_;

  //IndixCalculatorSection
  const std::size_t numModifications_;
  /** 
   * \brief Pointer to the first element in the points-(query-)Vector, 
   *  created for the scan,
   * This pointer is needed to recalculate the index of a specific element 
   * as we're passing pointers to the vector containing the queryIntervals.
   */
  const key_type * queryVectorPtrToFirstElement_;
  /** \brief Pointer to the first element in the interval-(data-)Vector, 
   *  created for the scan,
   *  this pointer is needed to recalculate the index of a specific element as 
   *  we're passing pointers to the vector containing the dataIntervals.
   */
  const key_type * dataVectorPtrToFirstElement_;

  //Randomizer section
  /** Random seed engine, has to be non-const as using the engine changes it. */
  std::mt19937 rSeedEngine_;
  typedef std::uniform_int_distribution<std::size_t> Distribution;

  const std::size_t offset_;
  const std::size_t cutoffSize_;
  /** Function pointer to a height calculator*/
  std::size_t (* const heightCalculator_)(const std::size_t);


 public:
 
  /** 
   * Constructor of the state we'll pass through most of the algorithm.
   *
   * \param limits The bounds in every dimension of the key, 
   *  they're needed to have values which are smaller/higher than all 
   *  possible values, see \ref HybridScanner
   * \param numMod The number of queries which are assigned to each input data: 
   *  As we want to calculate the true indices of the queries from their 
   *    pointers and several queries can belong to the same original object, 
   *    we have to use this.
   * \param queryVectorPtr Pointer to the first object in the queryVector, 
   *  needed to calculate the offset.
   * \param dataVectorPtr Pointer to the first object in the dataVector, 
   *  needed to calculate the offset.
   * \param cutoffSize Minimum amount of query and data objects to 
   *  handle in \ref HybridScanner, if one of them falls 
   *  below that value, switch to \ref OneWayScanner.
   * \param heightCalculator Function pointer to a heuristic which returns the 
   *  height to use in \ref getApproxMedian.
   */

  State(
      const key_type & limits, 
      const std::size_t numMod, 
      const key_type * queryVectorPtr,
      const key_type * dataVectorPtr,
      const std::size_t offset,
      const std::size_t cutoffSize = 250,
      std::size_t (*heightCalculator) (const std::size_t) = 
        &(SETA::State::defaultHeightCalculator_)
      ):
      limits_(limits), 
      numModifications_(numMod), 
      queryVectorPtrToFirstElement_(queryVectorPtr), 
      dataVectorPtrToFirstElement_(dataVectorPtr),
      offset_(offset),
      cutoffSize_(cutoffSize),
      heightCalculator_(heightCalculator)    
  {}

  inline 
  std::size_t calculate(bool isQueryVectorPtr, const key_type * objectPtr) const
  {
    if (isQueryVectorPtr)
    {
      return 
        offset_ + (objectPtr - this->queryVectorPtrToFirstElement_) / 
          this->numModifications_;
    }
    return (objectPtr - this->dataVectorPtrToFirstElement_);
  }

  inline std::size_t randInt(std::size_t lowerBound, std::size_t upperBound)
  {
    return Distribution(lowerBound, upperBound)(rSeedEngine_);
  }


  const key_type & getLimits() const { return limits_;}

  std::size_t heuristicHeight(const std::size_t n) const
  {
    return (*(this->heightCalculator_))(n);
  }
  std::size_t getCutoff() const{ return this->cutoffSize_; }
};




template <typename BoxType, std::size_t ...TIndices>
template <bool PointsContainQueries, std::size_t Dim>
struct SetA<BoxType, TIndices...>::
HybridScanner{

  /** 
   * \brief Find intersections by solving the 1-dimensional problems of 
   *  overlapping intervals.
   *
   * This is the main recursion function, it is a 
   *  divide & conquer algorithm, quite similar to quicksort.
   * By creating the nodes of a virtual segment tree in post-order, 
   *  we don't have to look at the whole tree at once (saving memory).
   * In 1 iteration, we can split the intervalset in 3 separate sets: 
   *   - Two of them represent the child-nodes of the virtual segment tree we're
   *      working on (we have to further seperate in the current dimension).
   *   - The third one is a set of solutions in the current dimension -> 
   *      we can solve the problem in the next dimension.
   *
   * \param[in] pointsPtrVector As we're looking at two subsets of intervals, 
   * this is the one representing the points.
   * \param[in] intervalsPtrVector These are the intervals, for this call.
   * \param[in] lowerBound As a recursion invariant, all points 
   *  represented by pointsPtrVector are inbetween the two 
   *  bounds (check recursion), which is why all intervals spanning across these
   *  bounds are intersecting with the points (in the current dimension).
   * \param[in] upperBound see lowerBound 
   * \param[in, out] state Contains a rng, can be used to track the 
   *  recursion and is able to calculate the indices.
   * \param[in, out] resultVector We pass the resultVector around to 
   *  add to it in OneWayScan
   * \note The resultVector is only modified in OneWayScanner so the 
   *  algorithm has to end there. We're switching to the
   *  OneWayScanner, which is more of a brute-force approach, when 
   *  either of the two cutoff values is met: 
   *   - Only one dimension is left
   *   - Either set is small enough that a match in one dimension will probably 
   *      lead to a match in the others.
   */ 
  static void scan(
    std::vector<const key_type *> & pointsPtrVector, //Points
    std::vector<const key_type *> & intervalsPtrVector,  //Intervals
    const typename std::tuple_element<Dim, key_type>::type::first_type & lowerBound,
    const typename std::tuple_element<Dim, key_type>::type::first_type & upperBound,
    State & state,
    ResultType & resultVector
    ) {

    typedef typename std::tuple_element<Dim, key_type>::type Key;
    typedef typename std::tuple_element<Dim, comp_type>::type Comp;
    Comp less;

    if (
      pointsPtrVector.size() == 0 || 
      intervalsPtrVector.size() == 0 || 
      !less(lowerBound, upperBound) 
    ) {
      return;
    }
    // switch into scanning mode if set sizes fall under the threshold
    if (
      pointsPtrVector.size() < state.getCutoff() || 
      intervalsPtrVector.size() < state.getCutoff() 
    ) {
      sortContainerHead<Dim>(pointsPtrVector);
      sortContainerHead<Dim>(intervalsPtrVector);
      OneWayScanner<PointsContainQueries, Dim>::
        scan(pointsPtrVector, intervalsPtrVector, state, resultVector);
      return;
    }
    // Set sizes are still above the threshold. We follow a divide and conquer
    // scheme: determine the median in the current dimension and use the
    // value to split the intervals and points for the next (recursive) call.
    std::size_t heuristicHeight = state.heuristicHeight(intervalsPtrVector.size());
    // FIXME: explain why you draw from the interval vector. This is the major
    // contribution of the implementation!!!
    typename Key::first_type median = 
      getApproxMedian<Dim>(intervalsPtrVector, heuristicHeight, state, less);
    
    std::vector<const key_type *> 
      intervalsMiddle, intervalsLeft, intervalsRight;

    typename std::vector<const key_type *>::const_iterator intVectorIt = 
      intervalsPtrVector.begin();

    while (intVectorIt != intervalsPtrVector.end()){
      Key interval = getKey<Dim>(*intVectorIt);
      const key_type * intPtr = *intVectorIt;
      ++intVectorIt;
      if (!less(lowerBound, interval.first)) { 
        if (!less(interval.second, upperBound)) { 
          intervalsMiddle.push_back(intPtr);
          continue;
        }
        intervalsLeft.push_back(intPtr);
        if (less(median, interval.second)) {
          intervalsRight.push_back(intPtr);
        }
      } else {
        if (!less(interval.first,median)){
          intervalsRight.push_back(intPtr);
          continue;
        }
        intervalsLeft.push_back(intPtr);
        if (less(median, interval.second))
          intervalsRight.push_back(intPtr);
      }
    }

    auto dimLimits = std::get<Dim+1>(state.getLimits());

    HybridScanner<PointsContainQueries, Dim+1>::
      scan(
        pointsPtrVector, 
        intervalsMiddle, 
        dimLimits.first, 
        dimLimits.second, 
        state, 
        resultVector
      );
    HybridScanner<!PointsContainQueries, Dim+1>::
      scan(
        intervalsMiddle, 
        pointsPtrVector, 
        dimLimits.first, 
        dimLimits.second, 
        state, 
        resultVector
      );

    //intervalsMiddle.swap(std::vector<const key_type *>());
    intervalsMiddle.clear();
    std::vector<const key_type *>().swap(intervalsMiddle);

    std::vector<const key_type *> pointsLeft, pointsRight;
    // all points which are to the left of the 
    // median have to be entered in pointsLeft.
    typename std::vector<const key_type *>::const_iterator pntVectorIt = 
      pointsPtrVector.begin();

    while (pntVectorIt != pointsPtrVector.end())
    {
      typename Key::first_type point = getHead<Dim>(*pntVectorIt);
      if (less(point, median)) pointsLeft.push_back(*pntVectorIt);
      else  pointsRight.push_back(*pntVectorIt);
      ++pntVectorIt;
    }
    HybridScanner<PointsContainQueries,Dim>::
      scan(pointsLeft, intervalsLeft, lowerBound, median, state, resultVector);
    intervalsLeft.clear();
    std::vector<const key_type *>().swap(intervalsLeft);
    pointsLeft.clear();
    std::vector<const key_type *>().swap(pointsLeft);
    HybridScanner<PointsContainQueries, Dim>::
      scan(pointsRight, intervalsRight, median, upperBound, state, resultVector);

  }

}; //end struct HybridScanner



/** 
 * This is a specialization of the HybridScanner when there's only 
 * one dimension left to compare in: just pass the sets to the OneWayScanner
 *
 * 
 */
template <typename BoxType, std::size_t ...TIndices>
template <bool PointsContainQueries>
struct SetA<BoxType, TIndices...>::
HybridScanner<PointsContainQueries, SetA<BoxType, TIndices...>::LASTDIM> {

  typedef SetA<BoxType, TIndices...> SETA;
  enum
  {
    LASTDIM = SETA::LASTDIM
  };
  typedef SETA::key_type key_type;
  /** standard scan*/
  inline static void scan(
      std::vector<const key_type *> & pointsPtrVector,
      std::vector<const key_type *> & intervalsPtrVector,
      const typename std::tuple_element<LASTDIM, key_type>::type::first_type & lowerBound,
      const typename std::tuple_element<LASTDIM, key_type>::type::first_type & upperBound,
      SETA::State & state,
      SETA::ResultType & resultVector
      )
  {
    SETA::sortContainerHead<LASTDIM>(pointsPtrVector);
    SETA::sortContainerHead<LASTDIM>(intervalsPtrVector);
    SETA::OneWayScanner<PointsContainQueries, LASTDIM>::
        scan(pointsPtrVector, intervalsPtrVector, state, resultVector);
  }
}; //end struct HybridScanner specialization



template <typename BoxType, std::size_t ...TIndices>
template <bool PointsContainQueries, std::size_t Dim>
struct SetA<BoxType, TIndices...>::
OneWayScanner{
  /**
   * \brief Pass through two sorted vectors and look for matches accordingly.
   * \param[in] pointsPtrVector As we're looking at two subsets of intervals, 
   * this is the one representing the points.
   * \param[in] intervalsPtrVector These are the intervals, for this call.
   * \param[in, out] state We need the state (containing 2 pointers) to 
   *  calculate the correct indices.
   * \param[in, out] resultVector Add our results. 
   *  \note As the OneWayScanner isn't necessarily be called for the 
   *  last dimension only, we have to use the IntersectionTester to 
   *  check for intersections in the remaining dimensions.
   */
  static void scan(
      const std::vector<const key_type * > & pointsPtrVector, 
      const std::vector<const key_type * > & intervalsPtrVector,
      State & state,
      ResultType & resultVector 
      ) {
    typedef typename std::tuple_element<Dim, key_type>::type::first_type Key;
    typedef typename std::tuple_element<Dim, comp_type>::type Comp; 
    typedef typename std::vector<const key_type * >::const_iterator CIT;
    typedef std::multiset<const key_type * , lessTail<Dim> > SortTailSet;
    SortTailSet intervalsPtrSet;
    typedef typename SortTailSet::const_iterator SIT;

    
    if (intervalsPtrVector.empty())
      return;

    Comp less;
    CIT pntVectorIt = pointsPtrVector.begin();
    CIT intVectorIt = intervalsPtrVector.begin();


    while (pntVectorIt != pointsPtrVector.end()){

      const key_type * pntPtr = *pntVectorIt;
      ++pntVectorIt; //don't look at the same point again!
      //const Key point = getHead<Dim>(pntPtr);
      key_type point = *pntPtr;
      Key lowerBound = std::get<Dim>(point).first; 
      std::get<Dim>(point).second = lowerBound;


      CIT oldIntVectorIt = intVectorIt; 
      while ( intVectorIt != intervalsPtrVector.end())
      {
        //if this is true, the lower endpoint of the intervals is greater than 
        //the query point - it is impossible for the query point to be inside.
        if (less(lowerBound,getHead<Dim>(*intVectorIt))) break;
        ++intVectorIt;
      }
      //add all intervals that weren't in the intervalSet yet whose lower end 
      //is not higher than the query point, these are the viable intervals.
      intervalsPtrSet.insert(oldIntVectorIt, intVectorIt);

      //return an iterator to the first object whose upper endpoint is 
      //greater than the point
      SIT activeSetIt = intervalsPtrSet.upper_bound(&point);
      //erase all intervals whose upper endpoints aren't greater than the point.
      intervalsPtrSet.erase(intervalsPtrSet.begin(), activeSetIt);
      //iterate through the data intervals till we find one whose 
      //lower endpoint is too big
      //fill up the intervalSet with all intervals before. 

      //erase all intervals from the active set whose upper endpoint are
      //lower than the query point, these can't encompass it.

      //add all intersections to the results
      SIT intersectionSetIt = intervalsPtrSet.begin();
      SIT intersectionSetEnd = intervalsPtrSet.end();

      for(; intersectionSetIt != intersectionSetEnd; ++intersectionSetIt) {
        if (IntersectionTester<Dim+1, LASTDIM+1>::
              test(pntPtr, *intersectionSetIt) ) {
            
          std::size_t edgeHead = state.calculate(PointsContainQueries, pntPtr);
          std::size_t edgeTail = state.calculate(!PointsContainQueries, *intersectionSetIt);
          resultVector[edgeHead].insert(resultVector[edgeHead].end(),static_cast<IntType>(edgeTail));
          resultVector[edgeTail].insert(resultVector[edgeTail].end(),static_cast<IntType>(edgeHead));
        }
      } //end add all intersections to the results.
    } //end while qContainerIt != pointsContainer.end() 
    // do loop for every query point.
  } //end scan

}; //end struct OneWayScanner


/**
 * \brief Check in \f$ O(Limit - Dim) \f$ if two intervals intersect by 
 *  comparing in all dimensions till it's either false
 *  or all tests evaluated to true.
 *   
 *
 * \param x Pointer to first interval.
 * \param y Pointer to second interval.
 * \tparam KeyDim The dimension of the _key_ we're primarily working on, 
 * this can differ from the actual dimension in value_type.
 * \tparam HeadDimension The dimension in the value_type we're working on, 
 *  needed to seperate the leftover dimensions.
 * \tparam TailDimensions... As we're looking for intersections in the 
 *  KeyDim (or HeadDim), we have to check for every intersection-pair
 *  if they're also matching in the other dimensions, the 
 *  parameter is mostly used for easier recursion.
 * 
 */


template <typename BoxType, std::size_t ...TIndices>
template <std::size_t Dim, std::size_t Limit>
struct SetA<BoxType, TIndices...>::
IntersectionTester {
  static bool test(const key_type * x, const key_type * y)
  {
    bool result;
    typedef typename std::tuple_element<Dim, key_type>::type Key;
    typedef typename std::tuple_element<Dim, comp_type>::type Comp;
    Key firstKey = getKey<Dim>(x);
    Key secondKey = getKey<Dim>(y);
    Comp less;
    if (less(firstKey.first, secondKey.first)){
      if (less(secondKey.first, firstKey.second)) result = true;
      else result = false;
    }
    else {
      if (less(firstKey.first, secondKey.second)) result = true;
      else result = false;
    }

    //bool short-circuiting saves us from comparing useless dimensions
    return result &&  IntersectionTester<Dim+1, Limit>::test(x, y); 

  }
};

template <typename BoxType, std::size_t ...TIndices>
template <std::size_t Limit>
struct SetA<BoxType, TIndices...>::
IntersectionTester<Limit, Limit> {
  typedef SetA<BoxType, TIndices...>::key_type key_type;
  static bool test(const key_type * x, const key_type * y){ return true; }
};



template <typename BoxType, std::size_t ... TIndices>
template <std::size_t Dim>
struct SetA<BoxType, TIndices...>::
lessHead{

  bool inline operator() (const key_type * x, const key_type * y) const {
    return 
        getCompareFunctor<Dim>()(
            getHead<Dim>(x),
            getHead<Dim>(y)
            );
  }
  bool inline operator() (const key_type & x, const key_type & y) const  {
    return 
        getCompareFunctor<Dim>()(
            getHead<Dim>(x),
            getHead<Dim>(y)
            );
  }
};


template <typename BoxType, std::size_t ... TIndices>
template <std::size_t Dim>
struct SetA<BoxType, TIndices...>::
lessTail{
  bool inline operator() (const key_type * x, const key_type * y) const {
    return 
        getCompareFunctor<Dim>()(
            getTail<Dim>(x),
            getTail<Dim>(y)
            );
  }

  bool inline operator() (const key_type & x, const key_type & y) const {
    return 
        getCompareFunctor<Dim>()(
            getTail<Dim>(x),
            getTail<Dim>(y)
            );
  }
};


/**
 * Sometimes we would like to see if the key we're looking at has the
 * correct values.
 * \tparam I Dimension the print should start with.
 * \tparam Number of dimensions the key_type possesses.
 */
template <typename BoxType, std::size_t ... TIndices>
template <std::size_t I, std::size_t N>
struct SetA<BoxType, TIndices...>::
KeyPrinter
{
  /** print key */
  static void 
      print(const key_type& key) {
        std::cout.setf(std::ios::fixed, std::ios::floatfield);
        std::cout << "Dim" << I <<": "<< 
          getHead<I>(key) << "-" << 
          getTail<I>(key) << '\n';
        KeyPrinter<I+1,N>::print(key);
      }
};

/**
 * Sometimes we would like to see if the key we're looking at has the
 * correct values, template specialization to terminate the recursion.
 */
template <typename BoxType, std::size_t ... TIndices>
template <std::size_t N>
struct SetA<BoxType, TIndices...>::
KeyPrinter<N,N>
{
  typedef SetA<BoxType, TIndices...>::key_type key_type;
  /** print key */
  static void 
      print(const key_type & key) {
        return;
      }
};

} //end namespace hybridtree;



#endif
