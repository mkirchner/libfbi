/* $Id: fbi-test.cpp 1 2010-10-30 01:14:03Z mkirchner $
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

#include <iostream>
#include <string>
#include <utility>
#include <list>
#include <functional>

//c++0x includes
//#include <tuple>
#include "unittest.hxx"
// hack for testing purposes
#define private public
#include <fbi/fbi.h>
#undef private
#include <fbi/traits.h>
//#include <fbi/tuplegenerator.h>
#include <fbi/connectedcomponents.h>
#include <boost/cstdint.hpp>

using namespace vigra;


// This is a ValueType for testing, the user has to create structs similar to
// this one: the typedefs dim_type is needed, same is true for getInterval
// which should be available for all active dimensions.

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
//BOOST_STATIC_ASSERT( (boost::is_same<testtype, correcttype>::type::value));



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

#define MAKEPAIR(z,n,text) std::make_pair BOOST_PP_LPAREN() BOOST_PP_CAT(text,BOOST_PP_MUL(2,n)) BOOST_PP_COMMA() BOOST_PP_CAT(text, BOOST_PP_INC(BOOST_PP_MUL(2,n))) BOOST_PP_RPAREN()

template < BOOST_PP_ENUM_BINARY_PARAMS(MAX_DIMENSIONS, typename T, = boost::mpl::void_ BOOST_PP_INTERCEPT) >
struct ValueType {

typedef typename fbi::mpl::TraitsGenerator<BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, T)>::key_type key_type;
typedef typename fbi::mpl::TraitsGenerator<BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, T)>::single_key_type single_key_type;

key_type key_;

#define BOOST_PP_LOCAL_MACRO(n) \
template <BOOST_PP_ENUM_PARAMS(BOOST_PP_MUL(2,n), typename P)> \
ValueType(BOOST_PP_ENUM_BINARY_PARAMS(BOOST_PP_MUL(2,n), P, p) ){\
  key_ = boost::tuples::make_tuple(BOOST_PP_ENUM(n, MAKEPAIR, p));\
}\

#define BOOST_PP_LOCAL_LIMITS (1, MAX_DIMENSIONS)
#include BOOST_PP_LOCAL_ITERATE()

ValueType(const key_type & newKey):key_(newKey){}
};




template <typename ValueType>
struct ValueTypeStandardAccessor
{
  typedef ValueType BoxType;
  template <size_t N>
  typename boost::tuples::element<N, typename BoxType::key_type >::type
  get (const BoxType & box) const {
     return boost::tuples::get<N>(box.key_);
  }

};

template <typename ValueType>
struct OffsetQueryAccessor
{

  typedef ValueType BoxType;
  typedef typename BoxType::single_key_type single_key_type;
  single_key_type offset_;

#define  BOOST_PP_LOCAL_MACRO(n) \
  template <BOOST_PP_ENUM_PARAMS(n, typename Offset)>\
  OffsetQueryAccessor(BOOST_PP_ENUM_BINARY_PARAMS(n, Offset, offset)){\
    offset_ = boost::tuples::make_tuple(BOOST_PP_ENUM_PARAMS(n, offset));\
    }\

#define BOOST_PP_LOCAL_LIMITS (1, MAX_DIMENSIONS)
#include BOOST_PP_LOCAL_ITERATE()


  template <size_t N>
  typename boost::tuples::element<N, typename BoxType::key_type >::type
  get (const BoxType & box) const {
    typename boost::tuples::element<N, typename BoxType::key_type>::type val = boost::tuples::get<N>(box.key_);
    typename boost::tuples::element<N, single_key_type>::type offset = boost::tuples::get<N>(offset_);
  
    return std::make_pair(val.first + offset, val.second + offset) ;
  }

};


namespace fbi{
/*
template <typename ... Dimensions>
struct Traits<ValueType<Dimensions...> >: public mpl::TraitsGenerator<Dimensions...>{};
*/

template < BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS,typename T) >
struct Traits<ValueType<BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, T) > >: 
  public mpl::TraitsGenerator<BOOST_PP_ENUM_PARAMS(MAX_DIMENSIONS, T) > {};

}

struct HybridSetATestSuite : vigra::test_suite {
  HybridSetATestSuite() : vigra::test_suite("HybridSetA") {
    add(testCase(&HybridSetATestSuite::testSetAType));
    add(testCase(&HybridSetATestSuite::testSetBType));
    add(testCase(&HybridSetATestSuite::testCountFunctor));
    
    add(testCase(&HybridSetATestSuite::testCreateKeyTypes));
    //add(testCase(&HybridSetATestSuite::testAccessor)); 
    // add(testCase(&HybridSetATestSuite::testPtrCreator));
    //add(testCase(&HybridSetATestSuite::testCompAccessor));
    //add(testCase(&HybridSetATestSuite::testSortFunctor));
    //add(testCase(&HybridSetATestSuite::testOneWayScan)); 

    add(testCase(&HybridSetATestSuite::testTwoWayScanIntersect));
    //add(testCase(&HybridSetATestSuite::testTwoWayScan));
    add(testCase(&HybridSetATestSuite::testHybridScanDual)); 

    add(testCase(&HybridSetATestSuite::testHybridScanBig));

    add(testCase(&HybridSetATestSuite::testHybridScanTop));
    add(testCase(&HybridSetATestSuite::testHybridScanBottom));

    add(testCase(&HybridSetATestSuite::testHybridScanLeft));
    add(testCase(&HybridSetATestSuite::testHybridScanRight));
    add(testCase(&HybridSetATestSuite::testHybridScanTopRightCorner));
    add(testCase(&HybridSetATestSuite::testHybridScanTopLeftCorner));
    add(testCase(&HybridSetATestSuite::testHybridScanBottomRightCorner));
    add(testCase(&HybridSetATestSuite::testHybridScanBottomLeftCorner));
    add(testCase(&HybridSetATestSuite::testHybridScanOverlap));
    add(testCase(&HybridSetATestSuite::testHybridScanRightSide));
    add(testCase(&HybridSetATestSuite::testHybridScanTopSide));
    add(testCase(&HybridSetATestSuite::testHybridScanLeftSide));
    add(testCase(&HybridSetATestSuite::testHybridScanBottomSide));
    add(testCase(&HybridSetATestSuite::testHybridScanFrontSide));
    add(testCase(&HybridSetATestSuite::testHybridScanBackSide));
    add(testCase(&HybridSetATestSuite::testHybridScanNoMatch));
    add(testCase(&HybridSetATestSuite::testHybridScanOnlyPoints));
    add(testCase(&HybridSetATestSuite::testHybridScanAllPointsOutside));
    add(testCase(&HybridSetATestSuite::testHybridScanFunctorVectors));
  
  }

  //typedef std::pair<int, std::less<int> > IntDimension;
  //typedef std::pair<double, std::less<double> > DoubleDimension;


  void testSetAType() {

    typedef ValueType<int, double, std::string> Map;
    typedef fbi::SetA<Map, 1,0 > TTT;
    BOOST_MPL_ASSERT_RELATION(TTT::NUMDIMS,==, 2);
    typedef boost::tuples::tuple<std::pair<double, double>, std::pair<int, int> > corr_key_type;
    BOOST_MPL_ASSERT_MSG((boost::is_same<TTT::key_type, corr_key_type>::value), TTT_HAS_THE_WRONG_KEYTYPE, (TTT::key_type));
    typedef boost::tuples::tuple<std::less<double>, std::less<int> > corr_comp_type;
    BOOST_MPL_ASSERT_MSG((boost::is_same<TTT::comp_type, corr_comp_type>::value), TTT_HAS_THE_WRONG_COMPTYPE, (TTT::comp_type));
  }
  
  void testSetBType() {

    typedef ValueType<double, short, int> Map;
    typedef fbi::SetA<Map, 0,2 > TTT;
    typedef ValueType<double, int, float, int> Map2;
    typedef TTT::SetB<Map2, 0,3 > QQQ;
    BOOST_MPL_ASSERT_RELATION(QQQ::QNUMDIMS,==, 2);
    typedef boost::tuples::tuple<std::pair<double, double>, std::pair<int, int> > corr_key_type;
    BOOST_MPL_ASSERT_MSG((boost::is_same<TTT::key_type, corr_key_type>::value), TTT_HAS_THE_WRONG_KEYTYPE, (TTT::key_type));
    BOOST_MPL_ASSERT_MSG((boost::is_same<QQQ::qkey_type, corr_key_type>::value), QQQ_HAS_THE_WRONG_KEYTYPE, (QQQ::qkey_type));
    typedef boost::tuples::tuple<std::less<double>, std::less<int> > corr_comp_type;
    BOOST_MPL_ASSERT_MSG((boost::is_same<QQQ::qcomp_type, corr_comp_type>::value), QQQ_HAS_THE_WRONG_COMPTYPE, (QQQ::qcomp_type));
  }
  

  void testCountFunctor() {
    std::vector<int> a(2),b(3),c(4),d(5),e(6);
    
    boost::tuples::tuple<int, std::vector<int>, int> test1 = 
      boost::tuples::make_tuple(0, a, 0);
    boost::tuples::tuple<int, double, std::vector<int> > test2 =
      boost::tuples::make_tuple(0, 0, b);
    boost::tuples::tuple<std::vector<int>, double, double> test3 =
      boost::tuples::make_tuple(c, 0, 0);
    boost::tuples::tuple<int, std::vector<int> > test4 = 
      boost::tuples::make_tuple(0, d);
    boost::tuples::tuple<std::vector<int> > test5 = 
      boost::tuples::make_tuple(e);
    boost::tuples::tuple<int, double > test6 = 
      boost::tuples::make_tuple(0, 0);
  
    if(fbi::mpl::FunctorChecker::count(test1) != 4)
    {
      failTest("FunctorChecker test 1 failed");
    }
    if(fbi::mpl::FunctorChecker::count(test2) != 5)
    {
      failTest("FunctorChecker test 2 failed");
    }
    if(fbi::mpl::FunctorChecker::count(test3) != 6)
    {
      failTest("FunctorChecker test 3 failed");
    }
    if(fbi::mpl::FunctorChecker::count(test4) != 6)
    {
      failTest("FunctorChecker test 4 failed");
    }
    if(fbi::mpl::FunctorChecker::count(test5) != 6)
    {
      failTest("FunctorChecker test 5 failed");
    }
    if(fbi::mpl::FunctorChecker::count(test6) != 2)
    {
      failTest("FunctorChecker test 6 failed");
    }
    
  
  }







  void testCreateKeyTypes(){

    typedef ValueType<int, double, std::string> Map;


    typedef fbi::SetA<Map, 2, 0, 1> TTT;
    std::vector<Map> testVector;
    testVector.push_back(Map(1, 2, 3.0,4.0,std::string("a"), std::string("b")));
    testVector.push_back(Map(5,6,7.0,8.0,std::string("c"), std::string("d")));
    testVector.push_back(Map(9,10,11.0,12.0,std::string("e"), std::string("f")));
    testVector.push_back(Map(13,14,15.0,16.0,std::string("g"), std::string("h")));

    std::vector<TTT::key_type> intervalVector = 
        TTT::KeyCreator<Map, 2,0,1>::getVector(testVector, boost::tuples::make_tuple(ValueTypeStandardAccessor<Map>()) );


    if(TTT::getKey<0>(intervalVector[0]) != std::make_pair(std::string("a"), std::string("b") ))
    {
      failTest("the string dimension in the first entry in your intervalVector isn't right");
    }

    if(TTT::getKey<1>(intervalVector[0]) != std::make_pair(1,2))
    {
      failTest("the integer dimension in the first entry in your intervalVector isn't right");
    }
    if(TTT::getKey<1>(intervalVector[1]) != std::make_pair(5,6))
    {
      failTest("the integer dimension in the second entry in your intervalVector isn't right");
    }

    if(TTT::getKey<2>(intervalVector[2]) != std::make_pair(11.0,12.0))
    {
      failTest("the double dimension in the third entry in your intervalVector isn't right");
    }
  }



  void testTwoWayScanIntersect(){
    //First check if the other Dimensions get checked correctly. 
    typedef ValueType<double, double, double> Map;
    typedef fbi::SetA<Map,0, 1> TTT;

    std::vector<Map> testVector; 
    testVector.push_back(Map( 1.0, 3.0, 4.0, 6.0  ,5.0,10.0));
    testVector.push_back(Map( 0.0, 2.0, 2.0, 5.0  ,4.0,8.0));
    testVector.push_back(Map( 3.0, 7.0, 3.0, 4.0  ,11.0,13.0));
    testVector.push_back(Map(-5.0, 0.1, 2.0, 5.0  ,7.0,10.0));

    std::vector<TTT::key_type> testIntervalVector = TTT::KeyCreator<Map, 0,1>::getVector(testVector, boost::tuples::make_tuple(ValueTypeStandardAccessor<Map>()));
    std::vector<const TTT::key_type * > testPtrVector= TTT::createPtrVector(testIntervalVector); 

    std::vector<bool> correctBoolVector;
    correctBoolVector.push_back(true);
    correctBoolVector.push_back(true);
    correctBoolVector.push_back(false);
    correctBoolVector.push_back(false);

    correctBoolVector.push_back(true);
    correctBoolVector.push_back(true);
    correctBoolVector.push_back(false);
    correctBoolVector.push_back(true);

    correctBoolVector.push_back(false);
    correctBoolVector.push_back(false);
    correctBoolVector.push_back(true);
    correctBoolVector.push_back(false);

    correctBoolVector.push_back(false);
    correctBoolVector.push_back(true);
    correctBoolVector.push_back(false);
    correctBoolVector.push_back(true);

    std::vector<bool> testBoolVector;
    for (size_t i = 0; i < testPtrVector.size(); ++i)
      for (size_t j = 0; j < testPtrVector.size(); ++j)
        testBoolVector.push_back(TTT::IntersectionTester<0,2>::test(testPtrVector[i], testPtrVector[j]));

    for (size_t i = 0; i < correctBoolVector.size(); ++i)
    {
      if( correctBoolVector[i] != testBoolVector[i])
      { std::cout << "wrong BoolVectorEntry: " << i << std::endl;
        failTest("IntersectionTest didn't give correct booleans");
      }
    }
  }

  void testHybridScanDual() {
    typedef ValueType<double, int, std::vector<double> > Map;
    typedef ValueType<double, std::vector<double>, int> QMap;
    typedef fbi::SetA<Map, 1,2> TTT;
    typedef TTT::SetB<QMap, 2,1> QQQ;
    typedef TTT::ResultType ResultType;
    typedef ValueTypeStandardAccessor<Map> StandardFunctor;
    typedef ValueTypeStandardAccessor<QMap> StandardFunctor2;
    std::vector<Map> testVector;
    std::vector<QMap> queryVector;

    ResultType hybridResults = QQQ::intersect(testVector, StandardFunctor(), queryVector, StandardFunctor2());

}

  void testHybridScanBig(){
    typedef ValueType<int, double> Map;
    typedef fbi::SetA<Map, 0,1> TTT;
    typedef TTT::ResultType ResultType;
    typedef ValueTypeStandardAccessor<Map> StandardFunctor;

    std::vector<Map> testVector;

    const size_t intDimension = 100;
    const size_t doubleDimension = 100;

    for(size_t i = 0; i < intDimension; ++i){
      for (size_t j = 0; j < doubleDimension; ++j){
        testVector.push_back(Map(2 * i, 2 * i + 1, 5.0 * j, 5.0 * j+1));
      }
    }

    ResultType correctResults(testVector.size());

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      if ( i % doubleDimension != doubleDimension - 1)
        correctResults[i].insert((boost::uint_fast32_t)i+1);
      if (i % doubleDimension != 0)
        correctResults[i].insert ((boost::uint_fast32_t)i-1);
      if (i >= doubleDimension)
        correctResults[i].insert ((boost::uint_fast32_t)(i - doubleDimension));
      if (i < (intDimension - 1) * doubleDimension)
        correctResults[i].insert ((boost::uint_fast32_t)(i + doubleDimension));
    }


    TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), OffsetQueryAccessor<Map>(0,5.5),OffsetQueryAccessor<Map>(2,0.0));

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      ResultType::value_type::const_iterator it1, it3;
      for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
           it1 != hybridResults[i].end() || it3 != correctResults[i].end();
           ++it1, ++it3) {
        if (*it1 != *it3 ) {
          std::cout << "correct: " << i << "-" << *it3 << std::endl;
          std::cout << "hybridResults " << std::endl;
          std::cout << i <<"-"<< *it1 << " " << i << std::endl;
          failTest("hybridScan gave a wrong result");
        }
      }
    }
  }
  
  void testHybridScanTop() {

    typedef ValueType<int, double> Map;
    typedef fbi::SetA<Map, 0,1> TTT;
    typedef TTT::ResultType ResultType;
    typedef ValueTypeStandardAccessor<Map> StandardFunctor;
    typedef OffsetQueryAccessor<Map> IntervalMover;

    std::vector<Map> testVector;

    const size_t intDimension=100;
    const size_t doubleDimension=100;

    for(size_t i = 0; i < intDimension; ++i){
      for (size_t j = 0; j < doubleDimension; ++j){
        testVector.push_back(Map(3 * i, 3 * i + 1, 5.0 * j, 5.0 * j+1));
      }
    }



    ResultType correctResults(testVector.size());

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      if ( i % doubleDimension != doubleDimension - 1)
        correctResults[i].insert((boost::uint_fast32_t)i+1);
      if (i % doubleDimension != 0)
        correctResults[i].insert ((boost::uint_fast32_t)i-1);
    }

    TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(0,4.1));

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      ResultType::value_type::const_iterator it1, it3;
      for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
           it1 != hybridResults[i].end() || it3 != correctResults[i].end();
           ++it1, ++it3) {
        if (*it1 != *it3 ) {
          std::cout << "correct: " << i << "-" << *it3 << std::endl;
          std::cout << "hybridResults " << std::endl;
          std::cout << i <<"-"<< *it1 << " " << i << std::endl;
          failTest("hybridScan gave a wrong result");
        }
      }
    }
  }


  void testHybridScanBottom(){

    typedef ValueType<int, double> Map;
    typedef fbi::SetA<Map, 0,1> TTT;
    typedef TTT::ResultType ResultType;
    typedef ValueTypeStandardAccessor<Map> StandardFunctor;
    typedef OffsetQueryAccessor<Map> IntervalMover;


    std::vector<Map> testVector;

    const size_t intDimension=100;
    const size_t doubleDimension=100;

    for(size_t i = 0; i < intDimension; ++i){
      for (size_t j = 0; j < doubleDimension; ++j){
        testVector.push_back(Map(3 * i, 3 * i + 1, 5.0 * j, 5.0 * j+1));
      }
    }

    ResultType correctResults(testVector.size());

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      if ( i % doubleDimension != doubleDimension - 1)
        correctResults[i].insert((boost::uint_fast32_t)i+1);
      if (i % doubleDimension != 0)
        correctResults[i].insert ((boost::uint_fast32_t)i-1);
    }


    TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(0,5.9));

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      ResultType::value_type::const_iterator it1, it3;
      for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
           it1 != hybridResults[i].end() || it3 != correctResults[i].end();
           ++it1, ++it3) {
        if (*it1 != *it3 ) {
          std::cout << "correct: " << i << "-" << *it3 << std::endl;
          std::cout << "hybridResults " << std::endl;
          std::cout << i <<"-"<< *it1 << " " << i << std::endl;
          failTest("hybridScan gave a wrong result");
        }
      }
    }
  }


  void testHybridScanRight(){

    typedef ValueType<double,int> Map;
    typedef fbi::SetA<Map, 0,1> TTT;
    typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

    std::vector<Map> testVector;

    const size_t intDimension=100;
    const size_t doubleDimension=100;

    for(size_t i = 0; i < intDimension; ++i){
      for (size_t j = 0; j < doubleDimension; ++j){
        testVector.push_back(Map(3.0 * i, 3.0 * i + 1, 5 * j, 5 * j+1));
      }
    }



    ResultType correctResults(testVector.size());

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      if (i >= doubleDimension)
        correctResults[i].insert ((boost::uint_fast32_t)(i - doubleDimension));
      if (i < (intDimension - 1) * doubleDimension)
        correctResults[i].insert ((boost::uint_fast32_t)(i + doubleDimension));
    }


    TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(2.1,0));

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      ResultType::value_type::const_iterator it1, it3;
      for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
           it1 != hybridResults[i].end() || it3 != correctResults[i].end();
           ++it1, ++it3) {
        if (*it1 != *it3 ) {
          std::cout << "correct: " << i << "-" << *it3 << std::endl;
          std::cout << "hybridResults " << std::endl;
          std::cout << i <<"-"<< *it1 << " " << i << std::endl;
          failTest("hybridScan gave a wrong result");
        }
      }
    } 
  }
  void testHybridScanLeft(){

    typedef ValueType<int, double> Map;
    typedef fbi::SetA<Map, 0,1> TTT;
    typedef TTT::ResultType ResultType;
    typedef ValueTypeStandardAccessor<Map> StandardFunctor;
    typedef OffsetQueryAccessor<Map> IntervalMover;

    std::vector<Map> testVector;

    const size_t intDimension=100;
    const size_t doubleDimension=100;

    for(size_t i = 0; i < intDimension; ++i){
      for (size_t j = 0; j < doubleDimension; ++j){
        testVector.push_back(Map(10 * i, 10 * i + 3, 15.0 * j, 15.0 * j+3));
      }
    }




    ResultType correctResults(testVector.size());

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      if (i >= doubleDimension)
        correctResults[i].insert ((boost::uint_fast32_t)(i - doubleDimension));
      if (i < (intDimension - 1) * doubleDimension)
        correctResults[i].insert ((boost::uint_fast32_t)(i + doubleDimension));
    }

    TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(12,0.0));

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      ResultType::value_type::const_iterator it1, it3;
      for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
           it1 != hybridResults[i].end() || it3 != correctResults[i].end();
           ++it1, ++it3) {
        if (*it1 != *it3 ) {
          std::cout << "correct: " << i << "-" << *it3 << std::endl;
          std::cout << "hybridResults " << std::endl;
          std::cout << i <<"-"<< *it1 << " " << i << std::endl;
          failTest("hybridScan gave a wrong result");
        }
      }
    } 
  }
  void testHybridScanOverlap(){

    typedef ValueType<int, double> Map;
    typedef fbi::SetA<Map, 0,1> TTT;
    typedef TTT::ResultType ResultType;
    typedef ValueTypeStandardAccessor<Map> StandardFunctor;
    typedef OffsetQueryAccessor<Map> IntervalMover;

    std::vector<Map> testVector;

    const size_t intDimension=100;
    const size_t doubleDimension=100;

    for(size_t i = 0; i < intDimension; ++i){
      for (size_t j = 0; j < doubleDimension; ++j){
        testVector.push_back(Map(3 * i, 3 * i + 1, 5.0 * j, 5.0 * j+1));
      }
    }




    ResultType correctResults(testVector.size());

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      if (i >= doubleDimension)
        correctResults[i].insert ((boost::uint_fast32_t)(i - doubleDimension));
      if (i < (intDimension - 1) * doubleDimension)
        correctResults[i].insert ((boost::uint_fast32_t)(i + doubleDimension));
    }

    TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(3,0.0));

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      ResultType::value_type::const_iterator it1, it3;
      for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
           it1 != hybridResults[i].end() || it3 != correctResults[i].end();
           ++it1, ++it3) {
        if (*it1 != *it3 ) {
          std::cout << "correct: " << i << "-" << *it3 << std::endl;
          std::cout << "hybridResults " << std::endl;
          std::cout << i <<"-"<< *it1 << " " << i << std::endl;
          failTest("hybridScan gave a wrong result");
        }
      }
    }
  }
  void testHybridScanTopRightCorner(){

    typedef ValueType<int, double> Map;
    typedef fbi::SetA<Map, 0,1> TTT;
    typedef TTT::ResultType ResultType;
    typedef ValueTypeStandardAccessor<Map> StandardFunctor;
    typedef OffsetQueryAccessor<Map> IntervalMover;

    std::vector<Map> testVector;

    const size_t intDimension=100;
    const size_t doubleDimension=100;

    for(size_t i = 0; i < intDimension; ++i){
      for (size_t j = 0; j < doubleDimension; ++j){
        testVector.push_back(Map(10 * i, 10 * i + 3, 15.0 * j, 15.0 * j+3));
      }
    }




    ResultType correctResults(testVector.size());

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      if (i >= doubleDimension && i % doubleDimension != 0)
        correctResults[i].insert ((boost::uint_fast32_t)(i - doubleDimension - 1));
      if (i < (intDimension - 1) * doubleDimension && i % doubleDimension != doubleDimension-1)
        correctResults[i].insert ((boost::uint_fast32_t)(i + doubleDimension + 1));
    }
    TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(8,13.0));

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      ResultType::value_type::const_iterator it1, it3;
      for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
           it1 != hybridResults[i].end() || it3 != correctResults[i].end();
           ++it1, ++it3) {
        if (*it1 != *it3 ) {
          std::cout << "correct: " << i << "-" << *it3 << std::endl;
          std::cout << "hybridResults " << std::endl;
          std::cout << i <<"-"<< *it1 << " " << i << std::endl;
          failTest("hybridScan gave a wrong result");
        }
      }
    }  
  }
  void testHybridScanTopLeftCorner(){
    typedef ValueType<int, double> Map;
    typedef fbi::SetA<Map, 0,1> TTT;
    typedef TTT::ResultType ResultType;
    typedef ValueTypeStandardAccessor<Map> StandardFunctor;
    typedef OffsetQueryAccessor<Map> IntervalMover;

    std::vector<Map> testVector;

    const size_t intDimension=100;
    const size_t doubleDimension=100;

    for(size_t i = 0; i < intDimension; ++i){
      for (size_t j = 0; j < doubleDimension; ++j){
        testVector.push_back(Map(10 * i, 10 * i + 3, 15.0 * j, 15.0 * j+3));
      }
    }




    ResultType correctResults(testVector.size());

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      if (i >= doubleDimension && i % doubleDimension != 0)
        correctResults[i].insert ((boost::uint_fast32_t)(i - doubleDimension - 1));
      if (i < (intDimension - 1) * doubleDimension && i % doubleDimension != doubleDimension-1)
        correctResults[i].insert ((boost::uint_fast32_t)(i + doubleDimension + 1));
    }

    TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(12,13));

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      ResultType::value_type::const_iterator it1, it3;
      for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
           it1 != hybridResults[i].end() || it3 != correctResults[i].end();
           ++it1, ++it3) {
        if (*it1 != *it3 ) {
          std::cout << "correct: " << i << "-" << *it3 << std::endl;
          std::cout << "hybridResults " << std::endl;
          std::cout << i <<"-"<< *it1 << " " << i << std::endl;
          failTest("hybridScan gave a wrong result");
        }
      }
    } 
  }


  void testHybridScanBottomRightCorner(){

    typedef ValueType<int, double> Map;
    typedef fbi::SetA<Map, 0,1> TTT;
    typedef TTT::ResultType ResultType;
    typedef ValueTypeStandardAccessor<Map> StandardFunctor;
    typedef OffsetQueryAccessor<Map> IntervalMover;

    std::vector<Map> testVector;

    const size_t intDimension=100;
    const size_t doubleDimension=100;

    for(size_t i = 0; i < intDimension; ++i){
      for (size_t j = 0; j < doubleDimension; ++j){
        testVector.push_back(Map(10 * i, 10 * i + 3, 15.0 * j, 15.0 * j+3));
      }
    }




    ResultType correctResults(testVector.size());

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      if (i >= doubleDimension && i % doubleDimension != 0)
        correctResults[i].insert ((boost::uint_fast32_t)(i - doubleDimension - 1));
      if (i < (intDimension - 1) * doubleDimension && i % doubleDimension != doubleDimension-1)
        correctResults[i].insert ((boost::uint_fast32_t)(i + doubleDimension + 1));
    }

    TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(8, 17.0));

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      ResultType::value_type::const_iterator it1, it3;
      for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
           it1 != hybridResults[i].end() || it3 != correctResults[i].end();
           ++it1, ++it3) {
        if (*it1 != *it3 ) {
          std::cout << "correct: " << i << "-" << *it3 << std::endl;
          std::cout << "hybridResults " << std::endl;
          std::cout << i <<"-"<< *it1 << " " << i << std::endl;
          failTest("hybridScan gave a wrong result");
        }
      }
    }  }



    void testHybridScanBottomLeftCorner(){

      typedef ValueType<int, double> Map;
      typedef fbi::SetA<Map, 0,1> TTT;
      typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

      std::vector<Map> testVector;

      const size_t intDimension=100;
      const size_t doubleDimension=100;

      for(size_t i = 0; i < intDimension; ++i){
        for (size_t j = 0; j < doubleDimension; ++j){
          testVector.push_back(Map(10 * i, 10 * i + 3, 15.0 * j, 15.0 * j+3));
        }
      }




      ResultType correctResults(testVector.size());

      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        if (i >= doubleDimension && i % doubleDimension != 0)
          correctResults[i].insert ((boost::uint_fast32_t)(i - doubleDimension - 1));
        if (i < (intDimension - 1) * doubleDimension && i % doubleDimension != doubleDimension-1)
          correctResults[i].insert ((boost::uint_fast32_t)(i + doubleDimension + 1));
      }


      TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(12,17.0));

      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        ResultType::value_type::const_iterator it1, it3;
        for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
             it1 != hybridResults[i].end() || it3 != correctResults[i].end();
             ++it1, ++it3) {
          if (*it1 != *it3 ) {
            std::cout << "correct: " << i << "-" << *it3 << std::endl;
            std::cout << "hybridResults " << std::endl;
            std::cout << i <<"-"<< *it1 << " " << i << std::endl;
            failTest("hybridScan gave a wrong result");
          }
        }
      }  
    }

    void testHybridScanRightSide() {

      typedef ValueType<int, double,float> Map;
      typedef fbi::SetA<Map, 0,1,2> TTT;
      typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

      std::vector<Map> testVector;

      const size_t intDimension = 30;
      const size_t doubleDimension = 30;
      const size_t floatDimension = 30;

      for(size_t i = 0; i < intDimension; ++i) {
        for (size_t j = 0; j < doubleDimension; ++j) {
          for (size_t k = 0; k < floatDimension; ++k) {
            testVector.push_back(Map(3 * i, 3 * i + 2, 5.0 * j, 5.0 * j+1,4.0f * k, 4.0f * k+1 ));
          }
        }
      }



      ResultType correctResults(testVector.size());


      //for our test, we'll move the querybox by 1 in the x and y dimensions, we'll try to 
      //let the right side of the querybox meet the left side of a real box.
      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        if (i % floatDimension > 0)
          correctResults[i].insert ((boost::uint_fast32_t)(i - 1));
        if (i % floatDimension != floatDimension - 1)
          correctResults[i].insert ((boost::uint_fast32_t)(i + 1));
      }

      TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(-1,0.0,4.0f));

      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        ResultType::value_type::const_iterator it1, it3;
        for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
             it1 != hybridResults[i].end() || it3 != correctResults[i].end();
             ++it1, ++it3) {
          if (*it1 != *it3 ) {
            std::cout << "correct: " << i << "-" << *it3 << std::endl;
            std::cout << "hybridResults " << std::endl;
            std::cout << i <<"-"<< *it1 << " " << i << std::endl;
            failTest("hybridScan gave a wrong result");
          }
        }
      }  
    }

    void testHybridScanLeftSide() {

      typedef ValueType<int, double,float> Map;
      typedef fbi::SetA<Map, 0,1,2> TTT;
      typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

      std::vector<Map> testVector;

      const size_t intDimension = 30;
      const size_t doubleDimension = 30;
      const size_t floatDimension = 30;

      for(size_t i = 0; i < intDimension; ++i) {
        for (size_t j = 0; j < doubleDimension; ++j) {
          for (size_t k = 0; k < floatDimension; ++k) {
            testVector.push_back(Map(5* i,5*i + 3, 5.0 * j, 5.0 * j+1,4.0f * k, 4.0f * k+1 ));

          }
        }
      }



      ResultType correctResults(testVector.size());


      //for our test, we'll move the querybox by 1 in the x and y dimensions, we'll try to 
      //let the right side of the querybox meet the left side of a real box.
      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        if (i % floatDimension > 0)
          correctResults[i].insert ((boost::uint_fast32_t)(i - 1));
        if (i % floatDimension != floatDimension - 1)
          correctResults[i].insert ((boost::uint_fast32_t)(i + 1));
      }


      TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(2,0.0,4.0f));

      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        ResultType::value_type::const_iterator it1, it3;
        for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
             it1 != hybridResults[i].end() || it3 != correctResults[i].end();
             ++it1, ++it3) {
          if (*it1 != *it3 ) {
            std::cout << "correct: " << i << "-" << *it3 << std::endl;
            std::cout << "hybridResults " << std::endl;
            std::cout << i <<"-"<< *it1 << " " << i << std::endl;
            failTest("hybridScan gave a wrong result");
          }
        }
      }      
    }

    void testHybridScanTopSide(){

      typedef ValueType<int, double,float> Map;
      typedef fbi::SetA<Map, 0,1,2> TTT;
      typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

      std::vector<Map> testVector;

      const size_t intDimension = 30;
      const size_t doubleDimension = 30;
      const size_t floatDimension = 30;

      for(size_t i = 0; i < intDimension; ++i) {
        for (size_t j = 0; j < doubleDimension; ++j) {
          for (size_t k = 0; k < floatDimension; ++k) {
            testVector.push_back(Map(3 * i, 3 * i + 1, 5.0 * j, 5.0 * j + 3,4.0f * k, 4.0f * k+1 ));
          }
        }
      }



      ResultType correctResults(testVector.size());


      //for our test, we'll move the querybox by 1 in the x and y dimensions, we'll try to 
      //let the right side of the querybox meet the left side of a real box.
      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        if (i % floatDimension > 0)
          correctResults[i].insert ((boost::uint_fast32_t)(i - 1));
        if (i % floatDimension != floatDimension - 1)
          correctResults[i].insert ((boost::uint_fast32_t)(i + 1));
      }


      TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(0,-2.0,4.0f));

      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        ResultType::value_type::const_iterator it1, it3;
        for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
             it1 != hybridResults[i].end() || it3 != correctResults[i].end();
             ++it1, ++it3) {
          if (*it1 != *it3 ) {
            std::cout << "correct: " << i << "-" << *it3 << std::endl;
            std::cout << "hybridResults " << std::endl;
            std::cout << i <<"-"<< *it1 << " " << i << std::endl;
            failTest("hybridScan gave a wrong result");
          }
        }
      }      
    }

    void testHybridScanBottomSide(){

      typedef ValueType<int, double,float> Map;
      typedef fbi::SetA<Map, 0,1,2> TTT;
      typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

      std::vector<Map> testVector;

      const size_t intDimension = 30;
      const size_t doubleDimension = 30;
      const size_t floatDimension = 30;

      for(size_t i = 0; i < intDimension; ++i) {
        for (size_t j = 0; j < doubleDimension; ++j) {
          for (size_t k = 0; k < floatDimension; ++k) {
            testVector.push_back(Map(3 * i, 3 * i + 1, 5.0 * j, 5.0 * j + 3,4.0f * k, 4.0f * k+1 ));
          }
        }
      }



      ResultType correctResults(testVector.size());


      //for our test, we'll move the querybox by 1 in the x and y dimensions, we'll try to 
      //let the right side of the querybox meet the left side of a real box.
      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        if (i % floatDimension > 0)
          correctResults[i].insert ((boost::uint_fast32_t)(i - 1));
        if (i % floatDimension != floatDimension - 1)
          correctResults[i].insert ((boost::uint_fast32_t)(i + 1));
      }

      TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(0,2.0,4.0f));

      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        ResultType::value_type::const_iterator it1, it3;
        for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
             it1 != hybridResults[i].end() || it3 != correctResults[i].end();
             ++it1, ++it3) {
          if (*it1 != *it3 ) {
            std::cout << "correct: " << i << "-" << *it3 << std::endl;
            std::cout << "hybridResults " << std::endl;
            std::cout << i <<"-"<< *it1 << " " << i << std::endl;
            failTest("hybridScan gave a wrong result");
          }
        }
      }  
    }

    void testHybridScanFrontSide(){

      typedef ValueType<int, double,float> Map;
      typedef fbi::SetA<Map, 0,1,2> TTT;
      typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

      std::vector<Map> testVector;

      const size_t intDimension = 30;
      const size_t doubleDimension = 30;
      const size_t floatDimension = 30;

      for(size_t i = 0; i < intDimension; ++i) {
        for (size_t j = 0; j < doubleDimension; ++j) {
          for (size_t k = 0; k < floatDimension; ++k) {
            testVector.push_back(Map(3 * i, 3 * i + 1, 7.0 * j, 7.0 * j + 3, 5.0f*k, 5.0f*k+3 ));
          }
        }
      }



      ResultType correctResults(testVector.size());


      //for our test, we'll move the querybox by 1 in the x and y dimensions, we'll try to 
      //let the right side of the querybox meet the left side of a real box.
      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        if (i % floatDimension > 0)
          correctResults[i].insert ((boost::uint_fast32_t)(i - 1));
        if (i % floatDimension != floatDimension - 1)
          correctResults[i].insert ((boost::uint_fast32_t)(i + 1));
      }

      TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(0,0.,7.f));

      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        ResultType::value_type::const_iterator it1, it3;
        for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
             it1 != hybridResults[i].end() || it3 != correctResults[i].end();
             ++it1, ++it3) {
          if (*it1 != *it3 ) {
            std::cout << "correct: " << i << "-" << *it3 << std::endl;
            std::cout << "hybridResults " << std::endl;
            std::cout << i <<"-"<< *it1 << " " << i << std::endl;
            failTest("hybridScan gave a wrong result");
          }
        }
      }  
    }

    void testHybridScanBackSide() {

      typedef ValueType<int, double,float> Map;
      typedef fbi::SetA<Map, 0,1,2> TTT;
      typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

      std::vector<Map> testVector;

      const size_t intDimension = 30;
      const size_t doubleDimension = 30;
      const size_t floatDimension = 30;

      for(size_t i = 0; i < intDimension; ++i) {
        for (size_t j = 0; j < doubleDimension; ++j) {
          for (size_t k = 0; k < floatDimension; ++k) {
            testVector.push_back(Map(3 * i, 3 * i + 1, 7.0 * j, 7.0 * j + 2,5.0f * k, 5.0f * k+3 ));
          }
        }
      }



      ResultType correctResults(testVector.size());


      //for our test, we'll move the querybox by 1 in the x and y dimensions, we'll try to 
      //let the right side of the querybox meet the left side of a real box.
      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        if (i % floatDimension > 0)
          correctResults[i].insert ((boost::uint_fast32_t)(i - 1));
        if (i % floatDimension != floatDimension - 1)
          correctResults[i].insert ((boost::uint_fast32_t)(i + 1));
      }

      TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(0,0.0,3.0f));

      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        ResultType::value_type::const_iterator it1, it3;
        for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
             it1 != hybridResults[i].end() || it3 != correctResults[i].end();
             ++it1, ++it3) {
          if (*it1 != *it3 ) {
            std::cout << "correct: " << i << "-" << *it3 << std::endl;
            std::cout << "hybridResults " << std::endl;
            std::cout << i <<"-"<< *it1 << " " << i << std::endl;
            failTest("hybridScan gave a wrong result");
          }
        }
      }  
    }
    void testHybridScanNoMatch()
    {

      typedef ValueType<int, double,float> Map;
      typedef fbi::SetA<Map, 0,1,2> TTT;
      typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

      std::vector<Map> testVector;

      const size_t intDimension = 30;
      const size_t doubleDimension = 30;
      const size_t floatDimension = 30;

      for(size_t i = 0; i < intDimension; ++i) {
        for (size_t j = 0; j < doubleDimension; ++j) {
          for (size_t k = 0; k < floatDimension; ++k) {
            testVector.push_back(Map(4 * i, 4 * i + 1, 7.0 * j, 7.0 * j + 2,10.0f * k, 10.0f * k+3 ));
          }
        }
      }



      ResultType correctResults(testVector.size());
      TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(2,4.0,5.0f));

      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        ResultType::value_type::const_iterator it1, it3;
        for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
             it1 != hybridResults[i].end() || it3 != correctResults[i].end();
             ++it1, ++it3) {
          if (*it1 != *it3 ) {
            std::cout << "correct: " << i << "-" << *it3 << std::endl;
            std::cout << "hybridResults " << std::endl;
            std::cout << i <<"-"<< *it1 << " " << i << std::endl;
            failTest("hybridScan gave a wrong result");
          }
        }
      }  
    }

    void testHybridScanOnlyPoints()
    {

      typedef ValueType<int, double,float> Map;
      typedef fbi::SetA<Map, 0,1,2> TTT;
      typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

      std::vector<Map> testVector;

      const size_t intDimension = 100;
      const size_t doubleDimension = 100;
      const size_t floatDimension = 100;

      for(size_t i = 0; i < intDimension; ++i) {
        for (size_t j = 0; j < doubleDimension; ++j) {
          for (size_t k = 0; k < floatDimension; ++k) {
            testVector.push_back(Map(0,0,0.0,0.0,0.0f,0.0f));
          }
        }
      }



      ResultType correctResults(testVector.size());
      TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(0,0.0,0.0f));

      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        ResultType::value_type::const_iterator it1, it3;
        for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
             it1 != hybridResults[i].end() || it3 != correctResults[i].end();
             ++it1, ++it3) {
          if (*it1 != *it3 ) {
            std::cout << "correct: " << i << "-" << *it3 << std::endl;
            std::cout << "hybridResults " << std::endl;
            std::cout << i <<"-"<< *it1 << " " << i << std::endl;
            failTest("hybridScan gave a wrong result");
          }
        }
      }  
    }




    void testHybridScanAllPointsOutside()
    {
      typedef ValueType<int, double,float> Map;
      typedef fbi::SetA<Map, 0,1,2> TTT;
      typedef TTT::SetB<Map, 0,1,2> QQQ;
      typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

      std::vector<Map> testVector;

      const size_t intDimension = 100;
      const size_t doubleDimension = 100;
      const size_t floatDimension = 100;

      for(size_t i = 0; i < intDimension; ++i) {
        for (size_t j = 0; j < doubleDimension; ++j) {
          for (size_t k = 0; k < floatDimension; ++k) {
            testVector.push_back(Map(0,1,0.0,1.0,0.0f,1.0f));
          }
        }
      }

      ResultType correctResults(testVector.size());
      TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(1,0.0,0.0f));

      for (size_t i = 0; i < correctResults.size(); ++i)
      {
        ResultType::value_type::const_iterator it1, it3;
        for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
             it1 != hybridResults[i].end() || it3 != correctResults[i].end();
             ++it1, ++it3) {
          if (*it1 != *it3 ) {
            std::cout << "correct: " << i << "-" << *it3 << std::endl;
            std::cout << "hybridResults " << std::endl;
            std::cout << i <<"-"<< *it1 << " " << i << std::endl;
            failTest("hybridScan gave a wrong result");
          }
        }
      }  
    }


void testHybridScanFunctorVectors() 
{

    typedef ValueType<double, int> Map;
    typedef fbi::SetA<Map, 0,1> TTT;
    typedef TTT::ResultType ResultType;
      typedef ValueTypeStandardAccessor<Map> StandardFunctor;
      typedef OffsetQueryAccessor<Map> IntervalMover;

    std::vector<Map> testVector;

    const size_t intDimension=100;
    const size_t doubleDimension=100;

    for(size_t i = 0; i < intDimension; ++i){
      for (size_t j = 0; j < doubleDimension; ++j){
        testVector.push_back(Map(3.0 * i, 3.0 * i + 1, 5 * j, 5 * j+2));
      }
    }



    ResultType correctResults(testVector.size());

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      if (i % doubleDimension > 0)
        correctResults[i].insert ((boost::uint_fast32_t)(i - 1));
      if (i % doubleDimension < doubleDimension - 1)  
        correctResults[i].insert ((boost::uint_fast32_t)(i + 1));
      correctResults[i].insert((boost::uint_fast32_t)i);
      if (i >= doubleDimension)
        correctResults[i].insert ((boost::uint_fast32_t)(i - doubleDimension));
      if (i < (intDimension - 1) * doubleDimension)
        correctResults[i].insert ((boost::uint_fast32_t)(i + doubleDimension));
    }

    std::vector<IntervalMover> movers;
    movers.push_back(IntervalMover(0.0,0));
    movers.push_back(IntervalMover(2.1,0));


    TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), 
    movers, IntervalMover(0,4)
    );

    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      ResultType::value_type::const_iterator it1, it3;
      for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
           it1 != hybridResults[i].end() || it3 != correctResults[i].end();
           ++it1, ++it3) {
        if (*it1 != *it3 ) {
          std::cout << "correct: " << i << "-" << *it3 << std::endl;
          std::cout << "hybridResults " << std::endl;
          std::cout << i <<"-"<< *it1 << " " << std::endl;
         // failTest("hybridScan gave a wrong result");
        }
      }
    } 


}






}; //end HybridSetATestSuite


struct HybridSetAProfileSuite : vigra::test_suite {
  HybridSetAProfileSuite() : vigra::test_suite("HybridSetAProfiler")
  {
    add(testCase(&HybridSetAProfileSuite::testHybridScan));
  }

  void testHybridScan(){

    typedef ValueType<int, double,float> Map;
    typedef fbi::SetA<Map, 0,1,2> TTT;
    typedef TTT::ResultType ResultType;
    typedef ValueTypeStandardAccessor<Map> StandardFunctor;
    typedef OffsetQueryAccessor<Map> IntervalMover;

    std::vector<Map> testVector;

    const size_t intDimension = 100;
    const size_t doubleDimension = 100;
    const size_t floatDimension = 100;

    size_t l = 0;
    for(size_t i = 0; i < intDimension; ++i) {
      for (size_t j = 0; j < doubleDimension; ++j) {
        for (size_t k = 0; k < floatDimension; ++k) {
          testVector.push_back(Map(3 , 3+ 2, 7.0, 7.0 + 3.0,10.0f *l, 10.0f* l + 1  ));
          ++l;
        }
      }
    }



    ResultType correctResults(testVector.size());


    //for our test, we'll move the querybox by 1 in the x and y dimensions, we'll try to 
    //let the right side of the querybox meet the left side of a real box.
    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      correctResults[i].insert ((ResultType::value_type::value_type) i);
    }

    TTT::ResultType hybridResults = TTT::intersect(testVector,StandardFunctor(), IntervalMover(0,0.0,0.0f));
    // std::cout << correctResults.size() << " " << hybridResults.size() << std::endl;
    // std::cout << "finished" << std::endl;
    for (size_t i = 0; i < correctResults.size(); ++i)
    {
      ResultType::value_type::const_iterator it1, it3;
      for (it1 = hybridResults[i].begin(), it3 = correctResults[i].begin(); 
           it1 != hybridResults[i].end() || it3 != correctResults[i].end();
           ++it1, ++it3) {
        if (*it1 != *it3 ) {
          std::cout << "correct: " << i << "-" << *it3 << std::endl;
          std::cout << "hybridResults " << std::endl;
          std::cout << i <<"-"<< *it1 << " " << i << std::endl;
          failTest("hybridScan gave a wrong result");
        }
      }
    }  
  }
};

struct ConnectedComponentsTestSuite : vigra::test_suite {
  ConnectedComponentsTestSuite() : vigra::test_suite("ConnectedComponents")
  {
    add(testCase(&ConnectedComponentsTestSuite::testConnectedComponents));
  }

  void testConnectedComponents(){
    typedef int IntType;
    IntType nodes = 50;
    std::vector<std::vector<IntType> > adjacencyList(nodes);
    IntType clusterS = 5;

    for(IntType i = 0; i < nodes; ++i)
    {
      if (i % clusterS == clusterS-1) { 
        adjacencyList[i].push_back(i - (clusterS-1));
        adjacencyList[i - (clusterS-1)].push_back(i);
      }    
      else{
        adjacencyList[i].push_back(i + 1);
        adjacencyList[i+1].push_back(i);
      }
    }
    std::vector<IntType> labels;
    
    IntType nComponents = findConnectedComponents(adjacencyList, labels);
    if (nComponents != 10) 
      failTest("wrong Number of labels"); 
    IntType label = 0;
    for (IntType i = 0; i < nodes; ++i)
    {
      if (i % clusterS == 0) ++label;
      if (labels[i] != label) 
      {
        std::cout << i << ": " << labels[i] << std::endl;
        failTest("wrong labels");
      }
    }
  }

  
};









int main() {
  HybridSetATestSuite test;
  int success = test.run();
  std::cout << test.report() << std::endl;

  ConnectedComponentsTestSuite ccTest;
  int success1 = ccTest.run();
  std::cout << ccTest.report() << std::endl;

  HybridSetAProfileSuite profile;
  int success2 = profile.run();
  std::cout << profile.report() << std::endl;


  return success || success1 || success2;

  //return success || success1;
}

