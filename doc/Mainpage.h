/*!

\mainpage libfbi: a fast box intersection C++ library
\image html libfbi-small.png


\section sec_intro Introduction
\c libfbi computes intersections between two sets of boxes. This is an important
task in many data analysis algorithms that require e.g. neighborhood information
(in the form of an adjacency list) for all points in large datasets.

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

\section sec_install Installation
\c libfbi is a header-only library and does not require the installation of any
compiled libraries. 

\subsection sec_install_src Distribution Packages
\c libfbi distribution packages for are available for download from 
the \c libfbi project page at http://software.steenlab.org/libfbi.
To install, simply run the installer and ensure that the \c libfbi include
directory is in the include search path.

\subsection sec_install_src Building from Source
Because \c libfbi is a header-only library, there is no need to build the
library inorder to use it in software projects. 
The following illustrates how to build the examples, tests and
distribution packages. Bundling \c libfbi requires a working CMake build system
(available from http://cmake.org/) and CMake >= 2.6.

With cmake in the system path, the build process is
\verbatim
 tar xvzf libfbi-xxxxxxx.tar.gz
 mkdir libfbi-build
 cd libfbi-build
 cmake ../libfbi-xxxxxxx
 make
 make test
\endverbatim
And, optionally:
\verbatim
 make package
\endverbatim

*/

