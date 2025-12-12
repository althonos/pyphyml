cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

message(STATUS "Building ${SKBUILD_PROJECT_NAME} for ${CMAKE_SYSTEM_PROCESSOR}")

# --- Setup include folders ----------------------------------------------------

set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/include)

# --- Add scripts for Cython ---------------------------------------------------

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/scripts/cmake")

# include(FindSSE3)
# include(CythonExtension)
# include("src/scripts/cmake/FindAVX1.cmake")
# include("src/scripts/cmake/FindAVX2.cmake")
# include("src/scripts/cmake/FindAVX512.cmake")
# include("src/scripts/cmake/FindNEON.cmake")
# include("src/scripts/cmake/FindPOPCNT.cmake")
# include("src/scripts/cmake/FindSSE2.cmake")
# include("src/scripts/cmake/FindSSE4.cmake")
# include("src/scripts/cmake/CythonExtension.cmake")

# --- Compile code -------------------------------------------------------------

add_subdirectory("src")
