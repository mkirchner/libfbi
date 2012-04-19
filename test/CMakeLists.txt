# check if have valgrind
SET(HAVE_VALGRIND TRUE)
FIND_PROGRAM(VALGRIND_EXECUTABLE valgrind)
IF(VALGRIND_EXECUTABLE STREQUAL "VALGRIND_EXECUTABLE-NOTFOUND")
    MESSAGE(STATUS "WARNING: Could not find valgrind. Will NOT build memory tests.")
    SET(HAVE_VALGRIND FALSE)
ENDIF(VALGRIND_EXECUTABLE STREQUAL "VALGRIND_EXECUTABLE-NOTFOUND")

MACRO(ADD_LIBFBI_TEST name exe src)
    STRING(REGEX REPLACE "test_([^ ]+).*" "\\1" test "${exe}" )

    #build the test
    ADD_EXECUTABLE(${exe} ${src})
    
    #add test to global list of unit test
    ADD_TEST(${name} ${exe})
    MESSAGE(STATUS "Adding test for ${name}: ${exe}.")
    #Add target for the test
    STRING(REGEX REPLACE "test_([^ ]+).*" "unit_\\1" unittest_target "${exe}" )
    ADD_CUSTOM_TARGET(${unittest_target} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${exe})

    IF (HAVE_VALGRIND) 
        #add memory test using valgrind
        STRING(REGEX REPLACE "test_([^ ]+).*" "mem_\\1" memtest_name "${exe}" )
        LIST(APPEND memtest_names ${memtest_name})
        LIST(APPEND exe_names ${exe})
        #Add target for the memory test
        ADD_CUSTOM_TARGET(${memtest_name} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/memcheck.py ${CMAKE_CURRENT_BINARY_DIR}/${exe})
    ENDIF (HAVE_VALGRIND) 
    
    #Add target for the gcov coverage report for this test
    #IF(WITH_GCOV)
    #    STRING(REGEX REPLACE "test_([^ ]+).*" "gcov_\\1" gcov_target "${exe}" )
    #    ADD_CUSTOM_TARGET(${gcov_target} COMMAND /bin/bash ${LIBFBI_BINARY_DIR}/scripts/for_cmake/unit_coverage.sh ${test})
    #ENDIF(WITH_GCOV)
ENDMACRO(ADD_LIBFBI_TEST name exe src)


# include directories
INCLUDE_DIRECTORIES(
    ${LIBFBI_SOURCE_DIR}/test/include
    ${LIBFBI_SOURCE_DIR}/include
    ${LIBFBI_BINARY_DIR}/include
    ${LIBFBI_BINARY_DIR}/test/include
    ${CMAKE_CURRENT_SOURCE_DIR}
)

#### Sources
SET(SRCS_FBI fbi-test.cpp)

#### Tests
ADD_LIBFBI_TEST("fbi" test_fbi ${SRCS_FBI})

LIST(LENGTH memtest_names numtests)
IF(numtests GREATER 0)
    MATH(EXPR numtests ${numtests}-1)
    FOREACH(i RANGE 0 ${numtests})
        LIST(GET memtest_names ${i} memtest_name)
        LIST(GET exe_names ${i} exe)
        MESSAGE(STATUS "Adding memory test for ${memtest_name}: ${exe}.")
        ADD_TEST(${memtest_name}
            ${CMAKE_CURRENT_SOURCE_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/${exe} ${CMAKE_BINARY_DIR})
    ENDFOREACH(i RANGE 0 ${numtests}-1)
ENDIF(numtests GREATER 0)

#### Benchmark targets
ADD_CUSTOM_TARGET(benchmark-fbi
    ${LIBFBI_BINARY_DIR}/test/scripts/benchmark.sh fbi
    DEPENDS example-xic-construction
)
ADD_CUSTOM_TARGET(benchmark-kdtree
    ${LIBFBI_BINARY_DIR}/test/scripts/benchmark.sh kdtree
    DEPENDS kdtree-xic-construction
)


