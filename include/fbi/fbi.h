
#ifndef __LIBFBI_INCLUDE_FBI_WINDOWS_H__
#define __LIBFBI_INCLUDE_FBI_WINDOWS_H__

//C++
#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <utility>
#include <vector>
#include <iostream>

//boost
#include <boost/cstdint.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/count.hpp>
#include <boost/mpl/remove.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/back.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/transform.hpp>
#include<boost/mpl/placeholders.hpp>
#include <boost/mpl/apply.hpp>
#include <boost/mpl/lambda.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/quote.hpp>

//tree
#include <fbi/config.h>
#include <fbi/traits.h>
#include <fbi/tuplegenerator.h>

namespace fbi { 
  template <typename BoxType, BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(MAX_DIMENSIONS, int TIndex, -1) > 
    class SetA {
  private:
  SetA();
    /** Small self-referencing typedef */
    typedef SetA<BoxType, BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, TIndex)> SETA;
    /** Helper-struct, encapsulating most typedefs to extract necessary types */
    typedef fbi::mpl::indexFilter<BoxType, BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, TIndex)> typeHelper;
  
  /** A compile-time constant to mark recursion tails.*/
  enum {
    NUMDIMS = typeHelper::NUMDIMS
  };
  
  BOOST_MPL_ASSERT_RELATION(NUMDIMS, >, 0);

  typedef BoxType value_type;
  
   /**
   *  The key is represented by a pair of values for each dimension that 
   *  should be considered for intersection tests, 
   * these pairs are wrapped in a std::tuple for easy access via std::get.
   * \see \ref fbi::Traits
   * \note This key_type is a subset of the original key_type, 
   * depending on the indices we're using.
   */
  typedef typename typeHelper::key_type key_type;

 /**
   * As every dimension can have a different type, 
   * the appropriate comparison-operator
   * has to be known during compile-time. 
   * \see \ref fbi::Traits
   */
  typedef typename typeHelper::comp_type comp_type;
  



public:


  /** 
   * Every intersection between two elements will be represented by their 
   * corresponding index in the container holding BoxType, 
   * to save memory we decided to artificially 
   * limit the amount of elements which allows us
   * to use 32Bit integers instead of 64Bit.
   *
   */
  
  typedef boost::uint_fast32_t IntType;

  typedef std::vector<int>::size_type size_type;
  typedef std::vector<int>::difference_type diff_type;


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

  template <typename QueryType, BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(MAX_DIMENSIONS, int QIndex, -1) > 
  struct SetB;

  private:
/** 
    * \class KeyCreator
   * \brief Handle the creation of keys by using given functors. 
   * As we have to handle variable indices, both from SetA and SetB, we have to add 
   * another template class to use different int parameter packs.
   *
   */
  //template <int ...KeyCreatorIndices>
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
  template <bool PointsContainQueries, int Dim>
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
   * \tparam Dim where we should first look for an intersection, 
   * all following dimensions will be 
   * evaluated via comparisons between two objects.
   */
  template <bool PointsContainQueries, int Dim>
  struct OneWayScanner;

  /** 
   * \class lessHead
   * \brief This Functor compares the lower endpoints of 
   * two key_type objects in a given dimension.
   * \see std::pair<T>::operator<() 
   * \tparam Dim dimension to compare in.
   */ 
  template <int Dim>
  struct lessHead;

  /** 
   * \class lessTail
   * \brief This Functor compares the upper endpoints of 
   * two key_type objects in a given dimension.
   * \see std::pair<T>::operator<() 
   * \tparam Dim dimension to compare in.
   */

  template <int Dim>
  struct lessTail;

  /** 
   * \class IntersectionTester
   * \brief Check if two objects intersect in 
   * all dimensions in [Dim,Limit).
   *
   * If the amount of data is small enough, 
   * straightforward checking for two objects can be used.
   * Note that this should only be used after 
   * sufficient filtering of possible intersection candidates by the scanners.
   *
   */

  template <int Dim, int Limit>
  struct IntersectionTester;
 /** 
   * \class KeyPrinter 
   * \brief For Debug reasons, print all dimensions of a given key via cout.
   */
  template <int Dim, int Limit>
  struct KeyPrinter;
 
 
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
   * \param[in] dataContainer A container with a const_iterator, 
   * it has to hold the same type as the one the tree was instantiated with.
   * Both key sets will be extracted by using 
   * functors on every value in this container.
   * \param[in] ifunctor This has to be either a class with a public 
   * \verbatim get<Dim>(const BoxType & ) const; \endverbatim method or
   * a std::vector holding that kind of class, if multiple interval objects
   * per box should be made. 
   * The tree will create the keys by using the functor to extract 
   * the indexed dimensions from the BoxType.
   * \param[in] qfunctors Like ifunctor, 
   * these functors will each create a different key, they can also be of type
   * std::vector<functor>, for every object in each vector
   * one query will be created per box in dataContainer.
   * 
   * \return We'll return a std::vector<std::set<IntType> >, 
   * as parallel edges are possible but not desired for the end result.
   * An edge (pair of two boxes), consisting of 2 IntType values a and b
   * -which refer to the indices of these boxes in dataContainer-
   * is represented by b being in the set of tails belonging to the head a,
   * as there are |dataContainer| sets in the result.
   * Note that every edge is inserted twice as intersection is reflective and 
   * the result an undirected graph.
   *
   * \see \ref SetB::intersect
   * \callgraph
   */
 
#define BOOST_PP_LOCAL_MACRO(n) \
 template < \
 class BoxContainer,\
 typename IntervalFunctor,\
 BOOST_PP_ENUM_PARAMS(n, typename QueryFunctor)\
 >\
  static ResultType intersect(\
    const BoxContainer & dataContainer,\
    const IntervalFunctor & ifunctor,\
    BOOST_PP_ENUM_BINARY_PARAMS(n, const QueryFunctor, & qfunctor)\
    ) {\
      boost::tuples::tuple<BOOST_PP_ENUM_PARAMS(n,QueryFunctor)> qfunctors = boost::tuples::make_tuple(BOOST_PP_ENUM_PARAMS(n, qfunctor));\
      return SetB<BoxType,BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, TIndex)>::\
      wrappedIntersect(State::defaultCutoff, dataContainer, ifunctor, qfunctors);\
  }

#define BOOST_PP_LOCAL_LIMITS (1, MAX_QFUNCTORS)
#include BOOST_PP_LOCAL_ITERATE()


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
   * \param[in] cutoff If one of the containers has less elements than cutoff,
   * switch to brute-force algorithm.
   * \param[in] dataContainer A container with a const_iterator, 
   * it has to hold the same type as the one the tree was instantiated with.
   * Both key sets will be extracted by using 
   * functors on every value in this container.
   * \param[in] ifunctor This has to be either a class with a public 
   * \verbatim get<Dim>(const BoxType & ) const; \endverbatim method or
   * a std::vector holding that kind of class, if multiple interval objects
   * per box should be made. 
   * The tree will create the keys by using the functor to extract 
   * the indexed dimensions from the BoxType.
   * \param[in] qfunctors Like ifunctor, 
   * these functors will each create a different key, they can also be of type
   * std::vector<functor>, for every object in each vector
   * one query will be created per box in dataContainer.
   * 
   * \return We'll return a std::vector<std::set<IntType> >, 
   * as parallel edges are possible but not desired for the end result.
   * An edge (pair of two boxes), consisting of 2 IntType values a and b
   * -which refer to the indices of these boxes in dataContainer-
   * is represented by b being in the set of tails belonging to the head a
   * (and vice versa), as there are |dataContainer| sets in the result.
   * Note that every edge is inserted twice as intersection is reflective and 
   * the result an undirected graph.
   *
   * \see \ref SetB::intersect
   * \callgraph
   */

#define BOOST_PP_LOCAL_MACRO(n) \
 template < \
 class BoxContainer,\
 typename IntervalFunctor,\
 BOOST_PP_ENUM_PARAMS(n, typename QueryFunctor)\
 >\
  static ResultType thetaIntersect(\
    const size_type cutoff,\
    const BoxContainer & dataContainer,\
    const IntervalFunctor & ifunctor,\
    BOOST_PP_ENUM_BINARY_PARAMS(n, const QueryFunctor, & qfunctor)\
    ) {\
      boost::tuples::tuple<BOOST_PP_ENUM_PARAMS(n,QueryFunctor)> qfunctors = boost::tuples::make_tuple(BOOST_PP_ENUM_PARAMS(n, qfunctor));\
      return SetB<BoxType,BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, TIndex)>::\
      wrappedIntersect(cutoff, dataContainer, ifunctor, qfunctors);\
  }

#define BOOST_PP_LOCAL_LIMITS (1, MAX_QFUNCTORS)
#include BOOST_PP_LOCAL_ITERATE()




private:
 /**
   * Calculate an approximate median by using the median-of-three method
   * on a ternary tree of given height.
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
  template <int Dim>
  const 
  static 
  typename boost::tuples::element<Dim, key_type>::type::first_type
  getApproxMedian(
    const std::vector<const key_type *> & container,
    const int height, State & state, 
    const typename boost::tuples::element<Dim, comp_type>::type & less)
  {
    typedef typename boost::tuples::element<Dim, key_type>::type::first_type 
      ValType;
    
    if (height == 0) {
      int random = state.randInt(0, container.size()-1);
      int rand1 = state.randInt(0,1);
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
   * As our scanfunctions have to pass and copy containers of intervals, 
   * using a pointer based approach was deemed 
   * faster and safer as we won't modify the actual data.
   *
   * \see \ref KeyCreator::getVector()
   * \note Don't move/insert into the intervalContainer afterwards, as it will 
   * invalidate all pointers.
   */
  
  static const std::vector<const key_type *>
  createPtrVector (const std::vector<key_type> & container) {
    std::vector<const key_type *> ptrVector(container.size());
    typename std::vector<key_type>::const_iterator it = container.begin();
    int i = 0;
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

  /** Comfort function to get the interval in the correct dimension 
    \param[in] key The box object
  */
  template <int Dim>
  static inline typename boost::tuples::element<Dim,key_type>::type 
  getKey(const key_type & key) {
    return std::get<Dim>(key);
  }

  /** Comfort function to get the interval in the correct dimension 
    \param[in] key The box object
  */
  template <int Dim>
  static inline typename boost::tuples::element<Dim,key_type>::type 
  getKey(const key_type * key) {
    return std::get<Dim>(*key);
  }

  /** Comfort function to get the lower endpoint in the correct dimension 
    \param[in] key The box object
  */

  template <int Dim>
  static inline typename boost::tuples::element<Dim,key_type>::type::first_type 
  getHead(const key_type & key) {
    return getKey<Dim>(key).first;
  }

  /** Comfort function to get the lower endpoint in the correct dimension 
    \param[in] key The box object
  */
  template <int Dim>
  static inline typename boost::tuples::element<Dim,key_type>::type::first_type
  getHead(const key_type * key) {
    return getKey<Dim>(key).first;
  }

  /** Comfort function to get the upper endpoint in the correct dimension 
    \param[in] key The box object
  */
  template <int Dim>
  static inline typename boost::tuples::element<Dim,key_type>::type::first_type
  getTail(const key_type & key) {
    return getKey<Dim>(key).second;
  }



  /** Comfort function to get the upper endpoint in the correct dimension 
    \param[in] key The box object
  */
  template <int Dim>
  static inline typename boost::tuples::element<Dim,key_type>::type::first_type
  getTail(const key_type * key){
    return getKey<Dim>(key).second;
  }


  /** Comfort function to get the comparison functor in the correct dimension */
  template <int Dim>
  static inline const typename boost::tuples::element<Dim, comp_type >::type
  getCompareFunctor(){
    return typename boost::tuples::element<Dim,comp_type>::type();
  }


  /** Sort a container of pointers to key_type objects, compare their 
    * lower endpoints in the specified dimension 
    \param[in] container The container to sort.
  */ 
  template <int Dim>
  static inline void
  sortContainerHead(std::vector<const key_type * > & container){
    std::sort(container.begin(), container.end(), lessHead<Dim>());
  }

  /** Sort a container of pointers to key_type objects, compare their 
    * upper endpoints in the specified dimension 
    * \param[in] container The container to sort
    * 
    */ 
  template <int Dim>
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
   * \param[in] t1 first object to compare
   * \param[in] t2 second object to compare
   * \param[in] t3 third object to compare
   * \param[in] less comparison object
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
   * \param[in] numElements the height of the ternary tree shall be dependent
   * on the number of elements
   * \see \ref getApproxMedian() 
   */
  static int inline heuristicHeight(int numElements)
  {
    return (int)(log((double)numElements));
  }

}; //end class tree

template <typename BoxType, BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, int TIndex) > 
template <typename QBoxType, BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, int QIndex) > 
class SetA<BoxType,BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, TIndex)>::
SetB {

 private:


  /** Helper-struct, encapsulating most typedefs to extract necessary types */
  typedef fbi::mpl::indexFilter<QBoxType, BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, QIndex)> qTypeHelper;

  /** A compile-time constant to mark recursion tails.*/
  enum {
    QNUMDIMS = qTypeHelper::NUMDIMS
  };
  /** Ensure that the number of querydimensions for both sets are equal*/
  BOOST_MPL_ASSERT_RELATION(NUMDIMS, ==, QNUMDIMS);


  /** A comfort typedef to refer to the type of the input boxes*/
  typedef QBoxType qvalue_type;
  /**
   *  The key is represented by a pair of values for each dimension that 
   *  should be considered for intersection tests, 
   * these pairs are wrapped in a std::tuple for easy access via std::get.
   * \see \ref fbi::Traits
   * \note This key_type is a subset of the original key_type, 
   * depending on the indices we're using.
   */

typedef typename qTypeHelper::key_type qkey_type;
BOOST_MPL_ASSERT_MSG((boost::is_same<key_type, qkey_type>::value), KEY_TYPES_DONT_MATCH, (key_type, qkey_type));


 /**
   * As every dimension can have a different type, 
   * the appropriate comparison-operator
   * has to be known during compile-time. 
   * \see \ref fbi::Traits
   */
  typedef typename qTypeHelper::comp_type qcomp_type;      
    /** Ensure that the comparison operators are equal*/
  BOOST_MPL_ASSERT_MSG((boost::is_same<comp_type, qcomp_type>::value), KEY_TYPES_DONT_MATCH, (comp_type, qcomp_type));

 public:
  /**
   * \callgraph
   * \brief The public interface for the user to start the algorithm. 
   *
   * This function initializes the keys (consisting of intervals in 
   * several dimensions) of the intervalSet and the querySet, along with 
   * pointers to pass them to separation algorithms.
   * \ref State::calculate has to be used to 
   * recalculate the original indices from the pointers.
   *
   *
   * \param[in] dataContainer STL Container providing a 
   *  forward iterator and holding a value_type.
   * \param[in] ifunctor This has to be either a class with a public 
   * \verbatim get<Dim>(const BoxType & ) const; \endverbatim method or
   * a std::vector holding that kind of class, if multiple interval objects
   * per box should be made. 
   * The tree will create the keys by using the functor to extract 
   * the indexed dimensions from the BoxType. They will work on every object in
   * dataContainer.
   * \param[in] qdataContainer STL Container providing a 
   *  forward iterator and holding a qvalue_type.
   * \param[in] qfunctors Like ifunctor, 
   * these functors will each create a different key, they can also be of type
   * std::vector<functor>, for every object in each vector
   * one query will be created per object in qdataContainer.
   * (|qfunctors.size()| * |qdataContainer| for every qfunctors type)
   * \return We'll return a std::vector<std::set<IntType> >, 
   * as parallel edges are possible but not desired for the end result.
   * Note that in the bipartite case, the boxes in dataContainer will be indexed
   * from 0 upto |dataContainer|-1, whereas the ones in qdataContainer are 
   * identified by |dataContainer|..|dataContainer + qdataContainer|-1.
   * An edge (pair of two boxes), consisting of 2 IntType values a and b,
   * is represented by b being in the set of tails belonging to the head a
   * (and vice versa) as there are |dataContainer + qdataContainer|-1 sets in the 
   * result.
   * Note that every edge is inserted twice as intersection is reflective and 
   * the result an undirected graph.
   * 
   */
#define BOOST_PP_LOCAL_MACRO(n) \
 template < \
 class BoxContainer,\
 class QContainer,\
 typename IntervalFunctor,\
 BOOST_PP_ENUM_PARAMS(n, typename QueryFunctor)\
 >\
  static ResultType intersect(\
    const BoxContainer & dataContainer,\
    const IntervalFunctor & ifunctor,\
    const QContainer & qdataContainer,\
    BOOST_PP_ENUM_BINARY_PARAMS(n, const QueryFunctor, & qfunctor)\
    ) {\
      boost::tuples::tuple<BOOST_PP_ENUM_PARAMS(n,QueryFunctor)> qfunctors = boost::tuples::make_tuple(BOOST_PP_ENUM_PARAMS(n, qfunctor));\
      return SetB<BoxType,BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, TIndex)>::\
      wrappedIntersect(State::defaultCutoff, dataContainer, ifunctor, qfunctors);\
  }

#define BOOST_PP_LOCAL_LIMITS (1, MAX_QFUNCTORS)
#include BOOST_PP_LOCAL_ITERATE()

#define BOOST_PP_LOCAL_MACRO(n) \
 template < \
 class BoxContainer,\
 class QContainer,\
 typename IntervalFunctor,\
 BOOST_PP_ENUM_PARAMS(n, typename QueryFunctor)\
 >\
  static ResultType thetaIntersect(\
    const size_type cutoff,\
    const BoxContainer & dataContainer,\
    const IntervalFunctor & ifunctor,\
    const QContainer & qdataContainer,\
    BOOST_PP_ENUM_BINARY_PARAMS(n, const QueryFunctor, & qfunctor)\
    ) {\
      boost::tuples::tuple<BOOST_PP_ENUM_PARAMS(n,QueryFunctor)> qfunctors = boost::tuples::make_tuple(BOOST_PP_ENUM_PARAMS(n, qfunctor));\
      return SetB<BoxType,BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, TIndex)>::\
      wrappedIntersect(State::defaultCutoff, dataContainer, ifunctor, qfunctors);\
  }

#define BOOST_PP_LOCAL_LIMITS (1, MAX_QFUNCTORS)
#include BOOST_PP_LOCAL_ITERATE()



private:
  template <
  class BoxContainer,
        class QContainer,
        typename IntervalFunctor, 
        typename QueryFunctors
  > static
  ResultType wrappedIntersect(
      const size_type & cutoff,
      const BoxContainer & dataContainer, 
      const IntervalFunctor & ifunctor, 
      const QContainer & qdataContainer,
      const QueryFunctors& qfunctors
      ) {
        BOOST_MPL_ASSERT_RELATION(boost::tuples::length<qfunctors>::value, > 0);
    if (dataContainer.empty()) { return ResultType();}
    // Generate the set of query boxes. The BoxType is an arbitrary,
    // user-specified type, that does not necessarily have any notion of
    // dimensionality. This call converts the BoxType data into the 
    // K-dimenstional boxes for fast box intersection.
    const auto queryIntervalVector = KeyCreator<QIndices...>::
      getVector(qdataContainer, qfunctors);
    // Generate the set of data boxes. See above, just for the QueryBoxType.
    const auto dataIntervalVector = KeyCreator<TIndices...>::
      getVector(dataContainer, ifunctor);

    key_type limits = 
      make_tuple(
        std::get<TIndices>(Traits<value_type>::getLimits())
      ... );

    //const std::size_t numQueryFunctors = sizeof...(QueryFunctors); 
    
    const size_type numQueryFunctors = 
        mpl::FunctorChecker::count(qfunctors); 
    //if we're looking at two different sets, use different indices for the elements!
    const std::size_t offset = 
        (reinterpret_cast<const char* >(&(dataContainer)) == 
        reinterpret_cast<const char* >(&(qdataContainer))) ? 0 : dataContainer.size();
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

    auto dimLimits = boost::tuples::get<0>(state.getLimits()); 

    // Call the hybrid algorithm for stabbing queries in the interval vector.
    HybridScanner<true, NUMDIMS>::
      scan(
        pointsPtrVector, 
        intervalsPtrVector, 
        dimLimits.first, 
        dimLimits.second,state, 
        resultVector 
      );
    // Reverse the previous call: queries in the "point" vector.
    HybridScanner<false, NUMDIMS>::
      scan(
        intervalsPtrVector, 
        pointsPtrVector, 
        dimLimits.first, 
        dimLimits.second,state, 
        resultVector
      );
    return resultVector;
  }
//   /**
//   * \callgraph
//   * \brief The public interface for the user to start the algorithm. 
//   *
//   * This function initializes the keys (consisting of intervals in 
//   * several dimensions) of the intervalSet and the querySet, along with 
//   * pointers to pass them to separation algorithms.
//   * \ref State::calculate has to be used to 
//   * recalculate the original indices from the pointers.
//   *
//   * \param[in] cutoff The theta cutoff value for switching into OneWayScan
//   * \param[in] dataContainer STL Container providing a 
//   *  forward iterator and holding a value_type.
//   * \param[in] ifunctor This has to be either a class with a public 
//   * \verbatim get<Dim>(const BoxType & ) const; \endverbatim method or
//   * a std::vector holding that kind of class, if multiple interval objects
//   * per box should be made. 
//   * The tree will create the keys by using the functor to extract 
//   * the indexed dimensions from the BoxType. They will work on every object in
//   * dataContainer.
//   * \param[in] qdataContainer STL Container providing a 
//   *  forward iterator and holding a qvalue_type.
//   * \param[in] qfunctors Like ifunctor, 
//   * these functors will each create a different key, they can also be of type
//   * std::vector<functor>, for every object in each vector
//   * one query will be created per object in qdataContainer.
//   * (|qfunctors.size()| * |qdataContainer| for every qfunctors type)
//   * as parallel edges are possible but not desired for the end result.
//   * Note that in the bipartite case, the boxes in dataContainer will be indexed
//   * from 0 upto |dataContainer|-1, whereas the ones in qdataContainer are 
//   * identified by |dataContainer|..|dataContainer + qdataContainer|-1.
//   * An edge (pair of two boxes), consisting of 2 IntType values a and b,
//   * is represented by b being in the set of tails belonging to the head a
//   * (and vice versa) as there are |dataContainer + qdataContainer|-1 sets in the 
//   * result.
//   * Note that every edge is inserted twice as intersection is reflective and 
//   * the result an undirected graph.
//   */
//  template <
//  class BoxContainer,
//        typename = typename std::enable_if<std::is_same<typename BoxContainer::value_type, value_type>::value>::type,
//        class QContainer,
//        typename = typename std::enable_if<std::is_same<typename QContainer::value_type, qvalue_type>::value>::type,
//        typename IntervalFunctor, 
//        typename ... QueryFunctors
//  > static
//  ResultType thetaIntersect(
//      const size_t cutoff,
//      const BoxContainer & dataContainer, 
//      const IntervalFunctor & ifunctor, 
//      const QContainer & qdataContainer,
//      const QueryFunctors& ... qfunctors
//      ) {
//    static_assert( (sizeof...(QueryFunctors) > 0), 
//      "Need at least one query functor.");
//    if (dataContainer.empty()) { return ResultType();}
//    // Generate the set of query boxes. The BoxType is an arbitrary,
//    // user-specified type, that does not necessarily have any notion of
//    // dimensionality. This call converts the BoxType data into the 
//    // K-dimenstional boxes for fast box intersection.
//    const auto queryIntervalVector = KeyCreator<QIndices...>::
//      getVector(qdataContainer, qfunctors...);
//    // Generate the set of data boxes. See above, just for the QueryBoxType.
//    const auto dataIntervalVector = KeyCreator<TIndices...>::
//      getVector(dataContainer, ifunctor);
//
//    key_type limits = 
//      make_tuple(
//        std::get<TIndices>(Traits<value_type>::getLimits())
//      ... );
//
//    //const std::size_t numQueryFunctors = sizeof...(QueryFunctors); 
//    const std::size_t numQueryFunctors = 
//        mpl::FunctorChecker::count(qfunctors...); 
//    //if we're looking at two different sets, use different indices for the elements!
//    const std::size_t offset = 
//        (reinterpret_cast<const char* const>(&(dataContainer)) == 
//        reinterpret_cast<const char* const>(&(qdataContainer))) ? 0 : dataContainer.size();
//    State state(
//        limits,
//        numQueryFunctors,
//        &(queryIntervalVector[0]),
//        &(dataIntervalVector[0]),
//        offset,
//        cutoff
//        );
//
//    // Create a vector of pointers that reference the above query boxes. This
//    // allows us to work on pointers and save a bit of memory.
//    std::vector<const key_type *> pointsPtrVector = 
//      createPtrVector(queryIntervalVector);
//    std::vector<const key_type *> intervalsPtrVector = 
//      createPtrVector(dataIntervalVector);
//
//    ResultType resultVector(offset + qdataContainer.size());
//
//    auto dimLimits = std::get<0>(state.getLimits()); 
//
//    // Call the hybrid algorithm for stabbing queries in the interval vector.
//    HybridScanner<true, NUMDIMS>::
//      scan(
//        pointsPtrVector, 
//        intervalsPtrVector, 
//        dimLimits.first, 
//        dimLimits.second,state, 
//        resultVector 
//      );
//    // Reverse the previous call: queries in the "point" vector.
//    HybridScanner<false, NUMDIMS>::
//      scan(
//        intervalsPtrVector, 
//        pointsPtrVector, 
//        dimLimits.first, 
//        dimLimits.second,state, 
//        resultVector
//      );
//    return resultVector;
//  }
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
   * \see \ref createKey()
   */

  template <class Container, class ... Functors>
  static std::vector<key_type>
  getVector(const Container & container, const Functors& ...functors){
    /** If there is no functor, we can't extract data from the boxes.*/
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
   *  \see \ref createKey()
   */

  template <typename T>
  static void
  createKeys(typename std::vector<key_type>::iterator & intervalIt,
             const T & dataValue 
             ) {}

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
   * \see \ref createKey()
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
   *  Use functors on a given value_type object to 
   * create a key for every functor.
   *
   * \param intervalIt Iterator to a container which will 
   *  hold the different keys.
   * \param dataValue That's the box representing the class we're working on.
   * \param functor A vector of instances of Functor which provides a 
   *  std::pair<...> get<std::size_t>(value_type) member function
   * \param functors Other functors which will be used recursively.
   * \see \ref createKey()
   */
  template <typename T, typename Functor, typename ...Functors>
  static void 
  createKeys(
      typename std::vector<key_type>::iterator & intervalIt,
      const T & dataValue, 
      const std::vector<Functor> & functor, 
      const Functors & ...functors) 
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
   * \param[in] val The box object to create a key object with
   * \param[in] functor The functor which offers a templatized get method
   * to extract data out of val
   */
  template <typename T, typename Functor>
  static inline const
  key_type 
  createKey(const T& val, const Functor & functor) {
    return std::make_tuple(functor.template get<KeyCreatorIndices>(val)...);
  }

}; //end struct KeyCreator


template <typename BoxType, BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, int TIndex) > 
class SetA<BoxType,BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, TIndex)>::
State
{
 private: 
/** 
  * A reasonably good method to calculate the height of an approximate
  * median ternary tree
  * \param[in] n Number of elements to get the median from
  */
  static size_type defaultHeightCalculator_(const size_type n) {
    return 
      static_cast<size_type> (
        std::max( 0., -3.0 + 1.8 * log10( static_cast<double>( n ) ) )
      ); 
  }
  
  /** 
  * For every dimension, upper and lower bounds have to be
  * known, which will be larger or smaller than all
  * possible values.
  */
  const key_type limits_;

 /** 
  * As it is possible to create several keys per box,
  * we have to remember the number so index-calculations
  * can be done
  */
  const size_type numModifications_;
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
  /** We need a uniform distribution*/
  typedef std::uniform_int_distribution<size_type> Distribution;
/** As our second set of objects continues the
 * numbering scheme of the first, we have to add an offset
 * to the indices.
 */
  const size_type offset_;
 /** Also known as theta, the margin to switch to
  * brute-force 
  */
  const size_type cutoffSize_;
  /** Function pointer to a height calculator*/
  size_type (* const heightCalculator_)(const size_type);


 public:
 enum {
  defaultCutoff = 250
 };


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
   * \param offset Needed to correctly calculate the indices,
   * equal to sizeof(SetA)
   * \param cutoffSize Minimum amount of query and data objects to 
   *  handle in \ref HybridScanner, if one of them falls 
   *  below that value, switch to \ref OneWayScanner.
   * \param heightCalculator Function pointer to a heuristic which returns the 
   *  height to use in \ref getApproxMedian.
   */

  State(
      const key_type & limits, 
      const size_type numMod, 
      const key_type * queryVectorPtr,
      const key_type * dataVectorPtr,
      const size_type offset,
      const size_type cutoffSize = defaultCutoff,
      size_type (*heightCalculator) (const size_type) = 
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

 /** 
  * As we pass pointers around during our algorithm,
  * it is necessary to find the corresponding index
  * of the objects when writing pairs as results.
  * To do this, we first save the pointer to the
  * respective first object, using our knowledge of
  * the memory representation a simple pointer subtraction
  * yields the original index.
  * Note that for two different sets A and B, the boxes in B
  * get an offset of |A|.
  * \param[in] isQueryVectorPtr We have to keep track if
  * a pointer refers to a box in the first or second 
  * input set as that changes its index.
  * \param[in] objectPtr a pointer to a key, can be either
  * generated from the first or second type of boxes.
  */
  inline 
  size_type calculate(bool isQueryVectorPtr, const key_type * objectPtr) const
  {
    if (isQueryVectorPtr)
    {
      return 
        offset_ + (objectPtr - this->queryVectorPtrToFirstElement_) / 
          this->numModifications_;
    }
    return (objectPtr - this->dataVectorPtrToFirstElement_);
  }
/** 
 * In \ref getApproxMedian we have to pick a random object,
 * this returns a random integer between its two inputs.
 *
 * \param[in] lowerBound The returned integer will be
 * greater or equal to this
 * \param[in] upperBound The returned integer will be
 * lesser than this
 *
 *
 */
  inline size_type randInt(size_type lowerBound, size_type upperBound)
  {
    return Distribution(lowerBound, upperBound)(rSeedEngine_);
  }

  /** Getter */
  const key_type & getLimits() const { return limits_;}
 /** Return a good height for the ternary median tree
 * \param n Number of elements
 */
  size_type heuristicHeight(const size_type n) const
  {
    return (*(this->heightCalculator_))(n);
  }
  /** Getter*/
  size_type getCutoff() const{ return cutoffSize_; }
}; //end class State






} //end namespace fbi


#endif
