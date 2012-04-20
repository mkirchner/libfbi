ADD_SUBDIRECTORY(include)

INCLUDE_DIRECTORIES(
    ${LIBFBI_BINARY_DIR}/include
    ${LIBFBI_SOURCE_DIR}/examples/include
    ${LIBFBI_BINARY_DIR}/examples/include
    ${LIBFBI_SOURCE_DIR}/examples/include/spatial
    ${LIBFBI_BINARY_DIR}/examples/include/spatial
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${BOOST_INCLUDE_DIRS}
)

LINK_DIRECTORIES(
    ${BOOST_LIBRARY_DIRS}
)

ADD_EXECUTABLE(example-xic-construction example-xic-construction.cpp)
ADD_EXECUTABLE(example-isotope-patterns example-isotope-patterns.cpp)
ADD_EXECUTABLE(example-ms2-ms1-matching example-ms2-ms1-matching.cpp)
ADD_EXECUTABLE(simple-example simple-example.cpp)
ADD_EXECUTABLE(soccer-example soccer-example.cpp)
IF (HAS_KDTREE)
ADD_EXECUTABLE(kdtree-xic-construction kdtree-xic-construction.cpp)

TARGET_LINK_LIBRARIES(kdtree-xic-construction
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_DATE_TIME_LIBRARY}
)
ENDIF (HAS_KDTREE)
TARGET_LINK_LIBRARIES(example-ms2-ms1-matching
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_DATE_TIME_LIBRARY}
)

TARGET_LINK_LIBRARIES(example-xic-construction
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_DATE_TIME_LIBRARY}
)

TARGET_LINK_LIBRARIES(example-isotope-patterns
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_DATE_TIME_LIBRARY}
)

