#include "fbi/tuple.h"
#include "fbi/tuplegenerator.h"

#ifndef __LIBFBI_EXAMPLES_CENTROID_H__
#define __LIBFBI_EXAMPLES_CENTROID_H__
struct Centroid 
{
  double mz_;
  unsigned int sn_;
  double rt_;
  Centroid(const double& mz, const unsigned int  & sn, const double & rt) 
    : mz_(mz), sn_(sn), rt_(rt){}
};

namespace fbi {

  template<>
  struct Traits<Centroid> : mpl::TraitsGenerator<float, float> {};

} //end namespace fbi

struct CentroidBoxGenerator
{
  template <size_t N>
    typename fbi::tuple_element<N, 
             typename fbi::Traits<Centroid>::key_type>::type 
               get(const Centroid&) const;

  double mzOffset_;
  double mzWindowPpm_;
  double rtOffset_;
  double rtWindow_;
  float snWindow_;

  CentroidBoxGenerator(double mzWindowPpm, float snWindow)
    : mzOffset_(0.0), mzWindowPpm_(mzWindowPpm), snWindow_(snWindow)
  {}
  CentroidBoxGenerator(double mzWindowPpm, double rtWindow, unsigned int snWindow)
    : mzOffset_(0.0), mzWindowPpm_(mzWindowPpm),
    rtOffset_(0.0), rtWindow_(rtWindow), snWindow_(snWindow)
  {}

  CentroidBoxGenerator(double mzOffset, double mzWindowPpm, 
      double rtOffset, double rtWindow, unsigned int snWindow)
    : mzOffset_(mzOffset), mzWindowPpm_(mzWindowPpm),
    rtOffset_(rtOffset), rtWindow_(rtWindow), snWindow_(snWindow)
  {}


};


template <>
std::pair<float, float>  
CentroidBoxGenerator::get<0>(const Centroid & centroid) const 
{
  return std::make_pair(
      mzOffset_ + centroid.mz_* (1 - mzWindowPpm_ * 1E-6), 
      mzOffset_ + centroid.mz_* (1 + mzWindowPpm_ * 1E-6));
}


template <>
std::pair<float, float>  
CentroidBoxGenerator::get<1>(const Centroid & centroid) const 
{
  return std::make_pair(
      centroid.sn_ - snWindow_, 
      centroid.sn_ + snWindow_);
}





#endif
