PROJECT(LIBFBI)

# require at least cmake 2.6
cmake_minimum_required(VERSION 2.6)
# include macro definitions
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/macros/require_out_of_src_build.cmake)
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

#############################################################################
# versioning
#############################################################################
set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}")

#############################################################################
# build libraries
#############################################################################

# require out-of-source build
MACRO_REQUIRE_OUT_OF_SRC_BUILD("LIBFBI requires out-of-source builds.")

SET(CMAKE_CXX0X_FLAGS "-std=c++0x" )
# more warnings against everything
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion -Wall " )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX0X_FLAGS}" )

IF (NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
    ENDIF(NOT CMAKE_BUILD_TYPE)

TRY_COMPILE (HAS_UNIFORM_INT_DISTRIBUTION
    ${CMAKE_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/testrandom.cxx
)

CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fbi/config.h.cmake 
    ${LIBFBI_BINARY_DIR}/include/fbi/config.h
)

# add subdirs
OPTION(ENABLE_TESTING "Compile tests" ON)
OPTION(ENABLE_COVERAGE "Enable GCov coverage analysis (defines a 'coverage' target and enforces static build of ms++)" OFF)
OPTION(ENABLE_BENCHMARK "Do benchmarking and plotting" ON)
OPTION(ENABLE_EXAMPLES "Compile examples" ON)

IF (ENABLE_BENCHMARK OR ENABLE_EXAMPLES)
    # require boost
    SET(Boost_USE_STATIC_LIBS OFF)
    SET(Boost_USE_MULTITHREAD OFF)
    FIND_PACKAGE(Boost COMPONENTS program_options system REQUIRED)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
    LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
ENDIF (ENABLE_BENCHMARK OR ENABLE_EXAMPLES)

#############################################################################
# syncing docs
#############################################################################
ADD_CUSTOM_TARGET(syncdocs COMMAND ${LIBFBI_SOURCE_DIR}/scripts/sync-docs-gh-pages.sh ${LIBFBI_SOURCE_DIR} ${LIBFBI_BINARY_DIR})
ADD_DEPENDENCIES(syncdocs doc)

#############################################################################
# regression tests
#############################################################################
SET(UPDATE_TYPE git)
IF (ENABLE_TESTING)
  INCLUDE(CTest)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(test)
ENDIF (ENABLE_TESTING)

#############################################################################
# code coverage analysis
#############################################################################
IF (ENABLE_COVERAGE)
    SET (ENABLE_EXAMPLES ON)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    CONFIGURE_FILE(${LIBFBI_SOURCE_DIR}/cmake/scripts/coverage.sh.in
        ${LIBFBI_BINARY_DIR}/cmake/scripts/coverage.sh
        @ONLY IMMEDIATE
    )
    CONFIGURE_FILE(${LIBFBI_SOURCE_DIR}/cmake/scripts/unit_coverage.sh.in
        ${LIBFBI_BINARY_DIR}/cmake/scripts/unit_coverage.sh
        @ONLY IMMEDIATE
    )
    ADD_CUSTOM_TARGET(
        coverage COMMAND /bin/bash ${LIBFBI_BINARY_DIR}/cmake/scripts/coverage.sh
    )
ENDIF(ENABLE_COVERAGE)


#############################################################################
# time per number of nodes analysis
#############################################################################
IF (ENABLE_BENCHMARK)
    CONFIGURE_FILE(${LIBFBI_SOURCE_DIR}/test/testdata.zip
        ${LIBFBI_BINARY_DIR}/test/testdata.zip
        COPYONLY
    )
    CONFIGURE_FILE(${LIBFBI_SOURCE_DIR}/test/scripts/benchmark.sh.cmake
        ${LIBFBI_BINARY_DIR}/test/scripts/benchmark.sh
        @ONLY IMMEDIATE
    )
    CONFIGURE_FILE(${LIBFBI_SOURCE_DIR}/test/scripts/plotting.r.cmake
        ${LIBFBI_BINARY_DIR}/test/scripts/plotting.r
        COPYONLY
    )
    ADD_CUSTOM_TARGET(
        benchmark COMMAND /bin/bash ${LIBFBI_BINARY_DIR}/test/scripts/setupbenchmark.sh
    )
    ADD_DEPENDENCIES(benchmark centroid)
ENDIF(ENABLE_BENCHMARK)


############################################################################
# examples
############################################################################
IF (ENABLE_EXAMPLES)
    ADD_SUBDIRECTORY(examples)
ENDIF (ENABLE_EXAMPLES)



#############################################################################
# documentation
#############################################################################
# doxygen support
CONFIGURE_FILE(${LIBFBI_SOURCE_DIR}/doc/Doxyfile.cmake
    ${LIBFBI_BINARY_DIR}/doc/Doxyfile
    @ONLY IMMEDIATE
    )
FIND_PACKAGE(Doxygen)

# target for doxygen
ADD_CUSTOM_TARGET(doc COMMAND ${DOXYGEN} ${LIBFBI_BINARY_DIR}/doc/Doxyfile)

#############################################################################
# create a package config file
#############################################################################
# we don't need to bother about compiled libs - there are none.
SET (LIBFBI_INCLUDE ${LIBFBI_BINARY_DIR}/include)
CONFIGURE_FILE(
    ${LIBFBI_SOURCE_DIR}/LIBFBIConfig.cmake.in
    ${LIBFBI_BINARY_DIR}/LIBFBIConfig.cmake 
    @ONLY
)
CONFIGURE_FILE(
    ${LIBFBI_SOURCE_DIR}/LIBFBIUse.cmake.in
    ${LIBFBI_BINARY_DIR}/LIBFBIUse.cmake
    COPYONLY
)
INSTALL(FILES 
    ${LIBFBI_BINARY_DIR}/LIBFBIUse.cmake
    ${LIBFBI_BINARY_DIR}/LIBFBIConfig.cmake
    DESTINATION share/libfbi/cmake
    COMPONENT headers
)

#############################################################################
# installation
#############################################################################
# location, location, location
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
SET(CMAKE_INSTALL_PREFIX
"/usr/local" CACHE PATH "libfbi install prefix" FORCE
)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
# headers
INSTALL(DIRECTORY ${LIBFBI_SOURCE_DIR}/include/fbi
    DESTINATION include
    COMPONENT headers
    REGEX "/.git$" EXCLUDE
    PATTERN "config.h.cmake" EXCLUDE)

# config.h
INSTALL(FILES ${LIBFBI_BINARY_DIR}/include/fbi/config.h
    DESTINATION include/fbi
    COMPONENT headers)

# documentation
INSTALL(DIRECTORY ${LIBFBI_BINARY_DIR}/doc
    DESTINATION share/fbi
    COMPONENT documentation
    PATTERN "Doxyfile" EXCLUDE)

#############################################################################
# packaging
#############################################################################
# FIXME: make sure that the documentation is compiled and pdfs are 
#        generated prior to packaging!

set(CPACK_PACKAGE_NAME "libfbi")
set(CPACK_PACKAGE_VENDOR "steenlab.org")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Fast Box Intersection library")

# CPack version strings
SET(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})

SET(CPACK_RESOURCE_FILE_README "${LIBFBI_SOURCE_DIR}/README.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${LIBFBI_SOURCE_DIR}/LICENSE.txt")

set(CPACK_COMPONENTS_ALL libraries headers applications documentation)
set(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "Libraries")
set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "C++ Development Headers")
set(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "C++ Interface Documentation")

set(CPACK_COMPONENT_LIBRARIES_DESCRIPTION
    "Compiled libraries for fast box intersection.")
set(CPACK_COMPONENT_HEADERS_DESCRIPTION
    "C++ header files for libfbi development.")
set(CPACK_COMPONENT_DOCUMENTATION_DESCRIPTION
    "The C++ libfbi interface documentation (Doxygen).")

INCLUDE(CPack)

