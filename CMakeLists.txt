cmake_minimum_required (VERSION 3.6.3)
project (Pumas C)

# Build options
if (NOT CMAKE_BUILD_TYPE)
        set (CMAKE_BUILD_TYPE "Release" CACHE STRING
                 "Choose the type of build, options are: Debug Release
                  RelWithDebInfo MinSizeRel and Test." FORCE)
endif ()

if (NOT BUILD_SHARED_LIBS)
    set (BUILD_SHARED_LIBS TRUE CACHE BOOL "Build the library as shared")
endif ()

option (PUMAS_BUILD_TABULATE "Build the pumas-tabulate utility" OFF)

if ("${CMAKE_SYSTEM}" MATCHES "Linux")
        option (PUMAS_USE_GDB "Additional features for debugging with gdb" OFF)
endif ()

# Build and install rules for the PUMAS library.
add_library (pumas src/pumas.c include/pumas.h)
target_include_directories (pumas PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries (pumas m)

if (PUMAS_BUILD_TABULATE)
    target_compile_definitions (pumas PRIVATE "-D_BUILD_TABULATE")
endif ()

if (PUMAS_USE_GDB)
        target_compile_definitions (pumas PRIVATE "-DGDB_MODE")
endif ()

if (WIN32)
        target_compile_definitions (pumas PRIVATE "-D_CRT_SECURE_NO_WARNINGS")
    if (BUILD_SHARED_LIBS)
        target_compile_definitions (pumas
                PRIVATE "-DPUMAS_API=__declspec(dllexport)")
    endif ()
    install (TARGETS pumas DESTINATION .)
    install (FILES include/pumas.h DESTINATION .)
else ()
    install (TARGETS pumas DESTINATION lib)
    install (FILES include/pumas.h DESTINATION include)
endif ()

# Build and install rules for the tabulate utility, if enabled.
if (PUMAS_BUILD_TABULATE)
    add_executable (pumas-tabulate src/pumas-tabulate.c src/pumas-tabulate.h)
    target_link_libraries (pumas-tabulate pumas)
    if (WIN32)
        if (BUILD_SHARED_LIBS)
            target_compile_definitions (pumas-tabulate
                    PRIVATE "-DPUMAS_API=__declspec(dllimport)")
        endif ()
        install (TARGETS pumas-tabulate DESTINATION .)
    else ()
        install (TARGETS pumas-tabulate DESTINATION bin)
    endif ()
endif ()

# Build the test(s), if enabled
option (PUMAS_BUILD_TEST "Build the test suite" OFF)
if (PUMAS_BUILD_TEST)
    include ("tests/test-pumas.cmake")
endif ()
