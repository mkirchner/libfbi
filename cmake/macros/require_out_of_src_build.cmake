#
# Copyright (c) 2006, Alexander Neundorf, <neundorf@kde.org>
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#
macro(MACRO_REQUIRE_OUT_OF_SRC_BUILD errmsg)
   string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" building_in_src)
   if(building_in_src)
     message(SEND_ERROR "${errmsg}")
     message(FATAL_ERROR "Remove ${CMAKE_SOURCE_DIR}/CMakeCache.txt.")
   endif(building_in_src)
endmacro(MACRO_REQUIRE_OUT_OF_SRC_BUILD)

