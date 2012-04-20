/*!

\mainpage libfbi: fast box intersection in C++
\image html libfbi-small.png


\section sec_intro Introduction

\c libfbi is a header-only C++ template library that enables the efficient
solution of box intersection problems in an arbitrary number of dimensions. The
implementation is based on a publication by Zomorodian and Edelsbrunner (2000)
and makes heavy use of C++ metaprogramming and variadic template programming
techniques. Despite this complexity, the library provides a straightforward and
simple interface that allows easy integration into new and existing projects. 

libfbi has been developed to tackle feature extraction and data analysis tasks
for high-resolution liquid chromatorgraphy/mass spectrometry (LC/MS) data sets.
The library itself, however, is in no way limited to this application scenario:
it can easily be used to approach many kinds of multi-dimensional box
intersection and related problems (e.g. determining bounding volumes,
k-nearest-neighbor search, density estimation, correspondence estimation in
motion or pattern alignment tasks, etc.).

\section sec_license Licensing
\c libfbi is available under an MIT license. In brief, this means that
everybody, including academic institutions and (in particular) companies are
expliciltly allowed to make use of \c libfbi, as long as the original
disclaimers and copyright notices are left untouched. Companies are explicitly allowed to
incorporate \c libfbi into new or existing propietary software products and 
are not required to disclose any source code that uses the library (yes, you
can have the cake and eat it too).
Please check the license text for all details.

That said, all contributions to the library are much appreciated.

\section sec_citation Citation
If you use \c libfbi in a publication, you must cite the following paper:
\li Kirchner M*, Xu B*, Steen JAJ, Steen H (2011). libfbi: A C++ Implementation
for Fast Box Intersection and Application to Sparse Mass Spectrometry Data.
<i>Submitted</i>.

\c libfbi development has been carried out in the Steen & Steen Lab:
http://steenlab.org/

\section Software Requirements
Please note that \c libfbi makes heavy use of TR1 and C++0x features, in
particular tuples and variadic templates. Development was carried out under GNU
GCC 4.4, on an Ubuntu Linux system. Windows users currently need to make use of
the MinGW compiler as variadic template support under Microsoft Visual Studio
will only be available starting from Visual Studio 11.

\subsection Update
\c libfbi is now available for compilers without C++1x support via \c BOOST. CMake will, according to the capabilities of your compiler, provide the correct headers to use.
Note that there was a slight API change, include \c <fbi/tuple.h> instead of \verbatim <tuple>, along with \c fbi::tuple_element instead of \c std::tuple::element

\c libfbi itself is free of any third-party dependencies; the tests and
examples make use of the VIGRA \c unittest header file and depend on the BOOST
library. Memory tests require \c valgrind, coverage analysis depends on \c
gcov. 

\section sec_download Download
Tar archives with \c libfbi headers, a CMake config file (expecting the
installation in \c /usr/local) and the complete HTML documentation are
available from http://software.steenlab.org/libfbi

\section sec_install Installation
\c libfbi is a header-only library and does not require the installation of any
compiled libraries. To use libfbi, 
\li download a header file and documentation bundle, 
\li include the header files in your project and
\li adjust the compilation flags of your compiler to support C++0x functionality
(e.g. \c libfbi requires the \c -std=c++0x flag in gcc).

If you use the \c CMake build system, unpack the \c libfbi distribution to a
convenient location and include the following lines into your
\c CMakeLists.txt to use \c libfbi:
\verbatim
[...]
find_package(LIBFBI)
if (LIBFBI_FOUND)
    include (${LIBFBI_USE_FILE})
    add_executable(MyProgram MyProgram.cpp)
endif(LIBFBI_FOUND)
[...]
\endverbatim
If you include the \c share/libfbi/cmake directory of the \c libfbi
distribution in the project search path, you can simply build your project
using
\verbatim
cd <my project build dir>
cmake [optional cmake flags for your project] <my project source dir> 
\endverbatim
If you do not include the directory, you need to tell \c find_package where to
find the \c libfbi configuration information:
\verbatim
cd <my project build dir>
cmake -DLIBFBI_DIR=<install_dir>/share/libfbi/cmake \
    [optional cmake flags for your project] <my project source dir> 
\endverbatim

\subsection sec_install_src Building from Source
The following illustrates how to build the examples, tests and
distribution packages. Bundling \c libfbi requires a working CMake build system
(available from http://cmake.org/) and CMake >= 2.6.

With cmake in the system path, the build process is
\verbatim
 tar xvzf libfbi-xxxxxxx.tar.gz
 mkdir libfbi-build
 cd libfbi-build
 ccmake ../libfbi
 make
 make test
\endverbatim
And, optionally:
\verbatim
 make package
\endverbatim
*/

