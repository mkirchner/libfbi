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

#ifndef __LIBFBI_INCLUDE_FBI_CONFIG_H__
#define __LIBFBI_INCLUDE_FBI_CONFIG_H__

#cmakedefine HAS_UNIFORM_INT_DISTRIBUTION
#ifndef HAS_UNIFORM_INT_DISTRIBUTION
    #define uniform_int_distribution uniform_int
#endif

#cmakedefine HAS_VARIADIC_TEMPLATES
#ifdef HAS_VARIADIC_TEMPLATES
  #define __USE_VARIADIC_TEMPLATES__
#endif

#cmakedefine ENABLE_MULTITHREADING
#ifdef ENABLE_MULTITHREADING
  #define __LIBFBI_USE_MULTITHREADING__
#endif

#if _MSC_VER && !__INTEL_COMPILER
  #define __FBI_MSWORKAROUND__ 1
#else
  #define __FBI_MSWORKAROUND__ 0
#endif
#ifndef MAX_DIMENSIONS 
  #define MAX_DIMENSIONS 4
#endif
#ifndef MAX_QFUNCTORS
  #define MAX_QFUNCTORS 2
#endif

#endif 
