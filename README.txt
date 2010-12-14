================================================================================
                                 libfbi README
================================================================================

1. Summary
----------
libfbi is a header-only, C++ template library for fast box intersection queries
in an arbitrary number of dimensions.

2. Installation
---------------
libfbi comes as an installation package and installs its headers into standard
include directories (i.e. /usr/local/include on POSIX systems).
If that does not work for you, simply download the source distribution and copy
the header files from the libfbi include/ directory to a suitable location.
libfbi is a C++ header-only template library; there are no precompiled files.
Make sure that location is in your include search path and that your compiler
supports variadic templates (e.g. gcc >= 4.4). Additional details can
be found in the doc/ subdirectory and in the online documentation at
http://software.steenlab.org/libfbi.

3. Documentation
----------------
libfbi uses doxygen for source documentation and examples. The precompiled 
documentation is installed if you install an installation package. If you
downloaded a source distribution, build the 'doc' target to build the
documentation. Cutting edge online documentation is always available from
http://software.steenlab.org/libfbi.

