#ifndef __LIBFBI_INCLUDE_FBI_TUPLE_H__
#define __LIBFBI_INCLUDE_FBI_TUPLE_H__

#include <fbi/config.h>
#if defined(__USE_VARIADIC_TEMPLATES__)
#include <tuple>
namespace fbi{
#if defined(__USE_TEMPLATE_ALIASES__)
template<std::size_t __i, typename _Tp>
using tuple_element = std::tuple_element<__i, _Tp>;
#else //else __USE_TEMPLATE_ALIASES__
template <std::size_t __i, typename _Tp>
struct tuple_element {
  typedef typename std::tuple_element<__i, _Tp>::type type;
};
#endif //endif __USE_TEMPLATE_ALIASES__
} //end namespace fbi
#else 
#include <boost/tuple/tuple.hpp>
namespace fbi{
template <int N, class T>
struct tuple_element {
  typedef typename boost::tuples::element<N,T>::type type;
};
} //end namespace fbi

#endif 


#endif
