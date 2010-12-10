INCLUDE_DIRECTORIES(
    ${LIBFBI_SOURCE_DIR}/include
    ${LIBFBI_BINARY_DIR}/include
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

TARGET_LINK_LIBRARIES(example-ms2-ms1-matching
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
)

TARGET_LINK_LIBRARIES(example-xic-construction
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
)

TARGET_LINK_LIBRARIES(example-isotope-patterns
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
)

