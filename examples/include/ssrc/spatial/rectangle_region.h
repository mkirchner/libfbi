/*
 * Copyright 2003-2005 Daniel F. Savarese
 * Copyright 2006-2009 Savarese Software Research Corporation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.savarese.com/software/ApacheLicense-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __SSRC_SPATIAL_RECTANGLE_REGION_H
#define __SSRC_SPATIAL_RECTANGLE_REGION_H

#include <ssrc/libssrckdtree-packages.h>

__BEGIN_NS_SSRC_SPATIAL

template<typename Point,
         const unsigned int i = std::tr1::tuple_size<Point>::value - 1>
struct rectangle_region {
  Point lower, upper;

  //rectangle_region() = default;
  rectangle_region() : lower(), upper() { }

  rectangle_region(const Point & lower, const Point & upper) :
    lower(lower), upper(upper)
  { }

  // This is a helper function that is NOT part of the region oncept.
  static bool
  contains(const Point & point, const Point & lower, const Point & upper)
  {
    return ((std::tr1::get<i>(point) >= std::tr1::get<i>(lower) &&
             std::tr1::get<i>(point) <= std::tr1::get<i>(upper)) &&
            rectangle_region<Point, i-1>::contains(point, lower, upper));
  }

  bool contains(const Point & point) const {
    return rectangle_region::contains(point, lower, upper);
  }
};

template<typename Point>
struct rectangle_region<Point, 0> {
  Point lower, upper;

  //rectangle_region() = default;
  rectangle_region() : lower(), upper() { }

  rectangle_region(const Point & lower, const Point & upper) :
    lower(lower), upper(upper)
  { }

  // This is a helper function that is NOT part of the region oncept.
  static bool
  contains(const Point & point, const Point & lower, const Point & upper) {
    return !(std::tr1::get<0>(point) < std::tr1::get<0>(lower) ||
             std::tr1::get<0>(point) > std::tr1::get<0>(upper));
  }


  bool contains(const Point & point) const {
    return rectangle_region::contains(point, lower, upper);
  }
};

template<typename Point>
struct rectangle_region<Point, 1> {
  Point lower, upper;

  //rectangle_region() = default;
  rectangle_region() : lower(), upper() { }

  rectangle_region(const Point & lower, const Point & upper) :
    lower(lower), upper(upper)
  { }

  // This is a helper function that is NOT part of the region oncept.
  static bool
  contains(const Point & point, const Point & lower, const Point & upper) {
    return !(std::tr1::get<0>(point) < std::tr1::get<0>(lower) ||
             std::tr1::get<1>(point) < std::tr1::get<1>(lower) ||
             std::tr1::get<0>(point) > std::tr1::get<0>(upper) ||
             std::tr1::get<1>(point) > std::tr1::get<1>(upper));
  }

  bool contains(const Point & point) const {
    return rectangle_region::contains(point, lower, upper);
  }
};

template<typename Point>
struct rectangle_region<Point, 2> {
  Point lower, upper;

  //rectangle_region() = default;
  rectangle_region() : lower(), upper() { }

  rectangle_region(const Point & lower, const Point & upper) :
    lower(lower), upper(upper)
  { }

  // This is a helper function that is NOT part of the region oncept.
  static bool
  contains(const Point & point, const Point & lower, const Point & upper) {
    return !(std::tr1::get<0>(point) < std::tr1::get<0>(lower) ||
             std::tr1::get<1>(point) < std::tr1::get<1>(lower) ||
             std::tr1::get<2>(point) < std::tr1::get<2>(lower) ||
             std::tr1::get<0>(point) > std::tr1::get<0>(upper) ||
             std::tr1::get<1>(point) > std::tr1::get<1>(upper) ||
             std::tr1::get<2>(point) > std::tr1::get<2>(upper));
  }

  bool contains(const Point & point) const {
    return rectangle_region::contains(point, lower, upper);
  }
};

template<typename Point>
struct rectangle_region<Point, 3> {
  Point lower, upper;

  //rectangle_region() = default;
  rectangle_region() : lower(), upper() { }

  rectangle_region(const Point & lower, const Point & upper) :
    lower(lower), upper(upper)
  { }

  // This is a helper function that is NOT part of the region oncept.
  static bool
  contains(const Point & point, const Point & lower, const Point & upper) {
    return !(std::tr1::get<0>(point) < std::tr1::get<0>(lower) ||
             std::tr1::get<1>(point) < std::tr1::get<1>(lower) ||
             std::tr1::get<2>(point) < std::tr1::get<2>(lower) ||
             std::tr1::get<3>(point) < std::tr1::get<3>(lower) ||
             std::tr1::get<0>(point) > std::tr1::get<0>(upper) ||
             std::tr1::get<1>(point) > std::tr1::get<1>(upper) ||
             std::tr1::get<2>(point) > std::tr1::get<2>(upper) ||
             std::tr1::get<3>(point) > std::tr1::get<3>(upper));
  }

  bool contains(const Point & point) const {
    return rectangle_region::contains(point, lower, upper);
  }
};

__END_NS_SSRC_SPATIAL

#endif
