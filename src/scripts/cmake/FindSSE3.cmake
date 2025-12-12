#.rst:
# FindSSE3
# --------
#
# Finds SSE3 support
#
# This module can be used to detect SSE3 support in a C compiler.  If
# the compiler supports SSE3, the flags required to compile with
# SSE3 support are returned in variables for the different languages.
# The variables may be empty if the compiler does not need a special
# flag to support SSE3.
#
# The following variables are set:
#
# ::
#
#    SSE3_C_FLAGS - flags to add to the C compiler for SSE3 support
#    SSE3_FOUND - true if SSE3 is detected
#
#=============================================================================

set(_SSE3_REQUIRED_VARS)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${SSE3_FIND_QUIETLY})

# sample SSE3 source code to test
set(SSE3_C_TEST_SOURCE
"
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <pmmintrin.h>
#endif
int foo() {
    __m128 vOne = _mm_set1_ps(1.0);
    __m128 result = _mm_hsub_ps(vOne,vOne);
    return (int) _mm_cvtss_f32(result);
}
int main(void) { return foo(); }
")

# if these are set then do not try to find them again,
# by avoiding any try_compiles for the flags
if((DEFINED SSE3_C_FLAGS) OR (DEFINED HAVE_SSE3))
else()
  if(WIN32)
    set(SSE3_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts SSE3
      " "
      "/arch:SSE3")
  else()
    set(SSE3_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts SSE3
      " "
      #clang,GCC,GNU
      "-msse3"
    )
  endif()

  include(CheckCSourceCompiles)

  foreach(FLAG IN LISTS SSE3_C_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(HAVE_SSE3 CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try SSE3 C flag = [${FLAG}]")
    endif()
    check_c_source_compiles("${SSE3_C_TEST_SOURCE}" HAVE_SSE3)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(HAVE_SSE3)
      set(SSE3_C_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  unset(SSE3_C_FLAG_CANDIDATES)

  set(SSE3_C_FLAGS "${SSE3_C_FLAGS_INTERNAL}"
    CACHE STRING "C compiler flags for SSE3 intrinsics")
endif()

list(APPEND _SSE3_REQUIRED_VARS SSE3_C_FLAGS)

set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if(_SSE3_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(SSE3
                                    REQUIRED_VARS ${_SSE3_REQUIRED_VARS})

  mark_as_advanced(${_SSE3_REQUIRED_VARS})

  unset(_SSE3_REQUIRED_VARS)
else()
  message(SEND_ERROR "FindSSE3 requires C or CXX language to be enabled")
endif()
