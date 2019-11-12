cmake_minimum_required (VERSION 3.6.3)
project (PUMAS
    LANGUAGES C
    VERSION 0.15
)


# Build options
if (NOT CMAKE_BUILD_TYPE)
        set (CMAKE_BUILD_TYPE "Release" CACHE STRING
                 "Choose the type of build, options are: Debug Release
                  RelWithDebInfo MinSizeRel and Test." FORCE)
endif ()

if (NOT BUILD_SHARED_LIBS)
    set (BUILD_SHARED_LIBS TRUE CACHE BOOL "Build the library as shared")
endif ()

option (PUMAS_BUILD_EXAMPLES "Build the pumas example binaries" OFF)

if ("${CMAKE_SYSTEM}" MATCHES "Linux")
        option (PUMAS_USE_GDB "Additional features for debugging with gdb" OFF)
endif ()

if (WIN32)
    if (BUILD_SHARED_LIBS)
            set (PUMAS_API "-DPUMAS_API=__declspec(dllexport)")
    endif ()
    set (PUMAS_BIN ".")
else ()
    set (PUMAS_BIN "bin")
endif ()


# Build and install rules for the PUMAS library
add_library (pumas src/pumas.c include/pumas.h)
set_target_properties (pumas PROPERTIES VERSION ${PUMAS_VERSION})
target_include_directories (pumas PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_compile_definitions (pumas PRIVATE ${PUMAS_API} -DPUMAS_VERSION=${PUMAS_VERSION})
target_link_libraries (pumas m)

if (PUMAS_USE_GDB)
        target_compile_definitions (pumas PRIVATE "-DGDB_MODE")
endif ()

if (WIN32)
    target_compile_definitions (pumas PRIVATE "-D_CRT_SECURE_NO_WARNINGS")
    install (TARGETS pumas DESTINATION .)
    install (FILES include/pumas.h DESTINATION .)
else ()
    install (TARGETS pumas DESTINATION lib)
    install (FILES include/pumas.h DESTINATION include)
endif ()


# Build and install rules for the examples, if enabled
if (PUMAS_BUILD_EXAMPLES)
        macro (pumas_example __tag)
                set (__name "pumas-example-${__tag}")
                add_executable (${__name} "examples/${__tag}.c")
                target_compile_definitions (${__name} PRIVATE ${PUMAS_API})
                target_link_libraries (${__name} pumas)
                install (TARGETS ${__name} DESTINATION ${PUMAS_BIN} OPTIONAL)
        endmacro ()

        pumas_example (geometry)
        pumas_example (load)
        pumas_example (straight)
endif ()

# Build the test(s), if enabled
option (PUMAS_BUILD_TEST "Build the test suite" OFF)
if (PUMAS_BUILD_TEST)
    include ("tests/test-pumas.cmake")
endif ()
