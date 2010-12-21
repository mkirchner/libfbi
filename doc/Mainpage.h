/*!

\mainpage libfbi: a fast box intersection C++ library
\image html libfbi-small.png


\section sec_intro Introduction

\c libfbi is a header-only C++ template library that enables the efficient
solution of box intersection problems in an arbitrary number of dimensions. The
implementation is based on a publication by Zomorodian and Edelsbrunner (2000)
and makes heavy use of C++ metaprogramming and variadic template programming
techniques. Despite this complexity, the library provides a straightforward and
simple interface that allows easy integration into new and existing projects. 

libfbi has so far been applied to feature extraction and data analysis tasks
for high-resolution liquid chromatorgraphy/mass spectrometry (LC/MS) data sets.
The library itself, however, is in no way limited to this application scenario:
it can easily be used to tackle many kinds of multi-dimensional box
intersection and related problems (e.g. determining bounding volumes,
k-nearest-neighbor search, density estimation, correspondence estimation in
motion or pattern alignment tasks, etc.).

\section sec_license Licensing
\c libfbi is available under an MIT license. In essence, this means that
everybody, in particular including all academic institutions and companies are
expliciltly allowed to make use of \c libfbi, as long as the original
disclaimers and copyrights are left untouched. You are explicitly allowed to
incorporate \c libfbi into new or existing propietary software products and you
are not required to disclose any source code that uses the library or any
changes you make. For details, see the license
text available at the top of each header file and in LICENSE.txt.

That said, your bug reports, suggestions and bug fixes are appreciated.

\section sec_citation Citation
Details concerning the \c libfbi implementation are available from two principal
sources:
\li Kirchner M*, Xu B*, Steen JAJ, Steen H (2010). libfbi: A C++ Implementation
for Fast Box Intersection and Application to Sparse Mass Spectrometry Data.
<i>Submitted</i>.
\li Zomorodian A, Edelsbrunner H (2000). Fast Software for Box Intersections.
<i>Proceedings of the Sixteenth Annual Symposium on Computational Geometry,
129-138</i>.

\section Requirements
Please note that \c libfbi makes heavy use of TR1 and C++0x features, in
particular tuples and variadic templates. Development was carried out under GNU
GCC 4.4, on an Ubuntu Linux system. Windows users currently need to make use of
the MinGW compiler as variadic template support under Microsoft Visual Studio
will only be available starting from Visual Studio 11.

\c libfbi itself is free of any third-party dependencies; the tests and
examples make use of the VIGRA \c unittest header file and depend on the BOOST
library. Memory tests require \c valgrind, coverage analysis depends on \c
gcov. 

\section sec_install Installation
\c libfbi is a header-only library and does not require the installation of any
compiled libraries. To use libfbi, 
\li download a header file and documentation bundle, 
\li include the header files in your project and
\li adjust the compilation flags of your compiler (\c libfbi requires the 
\c -std=c++0x flag in gcc).

If you use the \c CMake build system, include the following lines into your
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

\subsection sec_install_src Building from Source
Because \c libfbi is a header-only library, there is no need to build any
library inorder to use it in software projects.
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

