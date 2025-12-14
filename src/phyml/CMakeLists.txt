# --- Configure ----------------------------------------------------------------

include(CheckIncludeFile)
include(CheckSymbolExists)
include(CheckLibraryExists)
include(CheckCCompilerFlag)
include(FindSSE3)

check_library_exists(m cos ""       HAVE_LIBM)

check_include_file("float.h"        HAVE_FLOAT_H)
check_include_file("inttypes.h"     HAVE_INTTYPES_H)
check_include_file("getopt.h"       HAVE_GETOPT_H)
check_include_file("stdint.h"       HAVE_STDINT_H)
check_include_file("stdio.h"        HAVE_STDIO_H)
check_include_file("stdlib.h"       HAVE_STDLIB_H)
check_include_file("strings.h"      HAVE_STRINGS_H)
check_include_file("string.h"       HAVE_STRING_H)
check_include_file("sys/stat.h"     HAVE_SYS_STAT_H)
check_include_file("sys/types.h"    HAVE_SYS_TYPES_H)
check_include_file("unistd.h"       HAVE_UNISTD_H)
check_include_file("vprintf.h"      HAVE_VPRINTF_H)

if(HAVE_LIBM)
    check_library_exists(m floor "" HAVE_FLOOR)
    check_library_exists(m pow   "" HAVE_POW)
    check_library_exists(m rint  "" HAVE_RINT)
    check_library_exists(m sqrt  "" HAVE_SQRT)
else()
    check_symbol_exists(floor "math.h" HAVE_FLOOR)
    check_symbol_exists(pow   "math.h" HAVE_POW)
    check_symbol_exists(rint  "math.h" HAVE_RINT)
    check_symbol_exists(sqrt  "math.h" HAVE_SQRT)
endif()

check_symbol_exists(strchr "string.h"   HAVE_STRCHR)
check_symbol_exists(strstr "string.h"   HAVE_STRSTR)

check_c_compiler_flag(-funroll-loops       HAVE_UNROLL_LOOPS)
check_c_compiler_flag(-fomit-frame-pointer HAVE_OMIT_FRAME_POINTER)
check_c_compiler_flag(-finline             HAVE_INLINE)
check_c_compiler_flag(-Winline             HAVE_WINLINE)

set(PACKAGE "phyml")
set(PACKAGE_BUGREPORT "stephane.guindon@lirmm.fr")
set(PACKAGE_NAME "PhyML")
set(PACKAGE_STRING "PhyML 3.3.20250429")
set(PACKAGE_TARNAME "phyml")
set(PACKAGE_URL "")
set(PACKAGE_VERSION "3.3.20250429")

# set(PHYML 1)
# set(PHYREX 0)
# set(PHYREXSIM 0)
# set(PHYTIME 0)
# set(RF 0)
# set(RWRAP 0)

set(STDC_HEADERS 1)
set(UNIX 1)
set(VERSION "3.3.20250429")

configure_file(
    "config.h.in"
    "config.h"
    NO_SOURCE_PERMISSIONS
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

# --- Copy/Patch files ---------------------------------------------------------

set(PHYML_SOURCES
    utilities.c 
    utilities.h
    optimiz.c 
    optimiz.h
    lk.c 
    lk.h
    bionj.c 
    bionj.h
    models.c 
    models.h
    free.c 
    free.h
    help.c 
    help.h
    simu.c 
    simu.h
    eigen.c 
    eigen.h
    pars.c 
    pars.h
    alrt.c 
    alrt.h
    interface.c 
    interface.h
    cl.c 
    cl.h
    spr.c 
    spr.h
    times.c 
    times.h
    draw.c 
    draw.h
    rates.c 
    rates.h
    mcmc.c 
    mcmc.h
    stats.c 
    stats.h
    mg.c 
    mg.h
    io.c 
    io.h
    make.c 
    make.h
    mixt.c
    mixt.h
    evolve.c 
    evolve.h
    init.c
    init.h
    nexus.c 
    nexus.h
    date.c 
    date.h
    tbe.c 
    tbe.h
    sse.c 
    sse.h
    avx.c 
    avx.h
    ancestral.c 
    ancestral.h
    xml.c 
    xml.h
    m4.c
    m4.h
    cv.c
    cv.h
    tiporder.h
    # tiporder.c
    invitee.h
    # invitee.c
    sse2neon.h
)

foreach(_file IN ITEMS ${PHYML_SOURCES})
    cmake_path(GET _file FILENAME _name)
    if(EXISTS ${PROJECT_SOURCE_DIR}/patches/${_name}.patch)
        add_custom_command(
            OUTPUT
                ${_file}
            COMMENT
                "Patching ${_file}"
            COMMAND
                ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py
                    --input ${PROJECT_SOURCE_DIR}/vendor/phyml/src/${_file}
                    --patch ${PROJECT_SOURCE_DIR}/patches/${_name}.patch
                    --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${PROJECT_SOURCE_DIR}/vendor/phyml/src/${_file}
                ${PROJECT_SOURCE_DIR}/patches/${_name}.patch
        )
    else()
        add_custom_command(
            OUTPUT
                ${_file}
            COMMENT
                "Copying ${_file}"
            COMMAND
                cmake -E copy ${PROJECT_SOURCE_DIR}/vendor/phyml/src/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${PROJECT_SOURCE_DIR}/vendor/phyml/src/${_file}
        )
    endif()
    set(PHYML_PATCHED_SOURCES ${PHYML_PATCHED_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endforeach()

# --- Build library ------------------------------------------------------------

add_library(phyml STATIC ${PHYML_PATCHED_SOURCES})
target_compile_definitions(phyml PUBLIC PHYML=1 DISABLE_CLI=1)
target_include_directories(phyml PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
set_property(TARGET phyml PROPERTY C_STANDARD 99)

if(HAVE_LIBM)
    target_link_libraries(phyml PUBLIC m)
endif()
if(SSE3_FOUND)
    target_compile_options(phyml PUBLIC ${SSE3_C_FLAGS})
    target_link_options(phyml PUBLIC ${SSE3_C_FLAGS})
endif()
if(HAVE_UNROLL_LOOPS)
    target_compile_options(phyml PRIVATE -funroll-loops)
endif()
if(HAVE_INLINE)
    target_compile_options(phyml PRIVATE -finline)
endif()
if(HAVE_WINLINE)
    target_compile_options(phyml PRIVATE -Winline)
endif()
