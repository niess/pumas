cmake_minimum_required (VERSION 2.8)
project (PUMAS C)

# Parse PUMAS library version
file(READ "include/pumas.h" _pumas_h)

string(REGEX MATCH "PUMAS_VERSION_MAJOR ([0-9]*)" _ ${_pumas_h})
set (PUMAS_VERSION_MAJOR ${CMAKE_MATCH_1})

string(REGEX MATCH "PUMAS_VERSION_MINOR ([0-9]*)" _ ${_pumas_h})
set (PUMAS_VERSION_MINOR ${CMAKE_MATCH_1})

string(REGEX MATCH "PUMAS_VERSION_PATCH ([0-9]*)" _ ${_pumas_h})
set (PUMAS_VERSION_PATCH ${CMAKE_MATCH_1})

set (PUMAS_VERSION_STRING ${PUMAS_VERSION_MAJOR}.${PUMAS_VERSION_MINOR}.${PUMAS_VERSION_PATCH})


# Build options
if (NOT CMAKE_BUILD_TYPE)
        set (CMAKE_BUILD_TYPE "Release" CACHE STRING
                 "Choose the type of build, options are: Debug Release
                  RelWithDebInfo and MinSizeRel." FORCE)
endif ()

if (NOT BUILD_SHARED_LIBS)
        set (BUILD_SHARED_LIBS TRUE CACHE BOOL "Build the library as shared")
endif ()

option (PUMAS_BUILD_EXAMPLES "Build the pumas example binaries" OFF)

if (UNIX)
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


# Disable extra Windows security checks when using msvc
if (MSVC)
        add_definitions("-D_CRT_SECURE_NO_WARNINGS")
endif ()


# Compile executables with a relative RPATH on Linux or OSX
if (PUMAS_BUILD_EXAMPLES)
        if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
                set (CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
        elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
                set (CMAKE_MACOSX_RPATH "ON")
                set (CMAKE_INSTALL_RPATH "@loader_path/../lib")
        endif ()
endif ()


if(WIN32)
	set(PUMAS_INCLUDE .)
	set(PUMAS_LIB .)
else()
	set(PUMAS_INCLUDE include)
	set(PUMAS_LIB lib)
endif()

# Build and install rules for the PUMAS library
add_library (pumas src/pumas.c include/pumas.h)
set_target_properties (pumas PROPERTIES
        VERSION ${PUMAS_VERSION_STRING}
        SOVERSION ${PUMAS_VERSION_MAJOR}
)
target_include_directories (pumas PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:${PUMAS_INCLUDE}>)
target_compile_definitions (pumas PRIVATE
        ${PUMAS_API}
        -DPUMAS_VERSION_MAJOR=${PUMAS_VERSION_MAJOR}
        -DPUMAS_VERSION_MINOR=${PUMAS_VERSION_MINOR}
        -DPUMAS_VERSION_PATCH=${PUMAS_VERSION_PATCH})
if (UNIX)
        target_link_libraries (pumas PUBLIC m)
endif ()

if (PUMAS_USE_GDB)
        target_compile_definitions (pumas PRIVATE "-DGDB_MODE")
endif ()

install (TARGETS pumas 
	EXPORT pumasTarget
	DESTINATION ${PUMAS_LIB})
install (FILES include/pumas.h DESTINATION ${PUMAS_INCLUDE})

install(EXPORT pumasTarget
	FILE pumasTarget.cmake
	DESTINATION "${PUMAS_LIB}/cmake/pumas")

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/pumasConfig.cmake"
	INSTALL_DESTINATION "${PUMAS_LIB}/cmake/pumas"
	NO_SET_AND_CHECK_MACRO
	NO_CHECK_REQUIRED_COMPONENTS_MACRO)

write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/pumasConfigVersion.cmake"
	VERSION "${PUMAS_VERSION_MAJOR}.${PUMAS_VERSION_MINOR}.${PUMAS_VERSION_PATCH}"
	COMPATIBILITY SameMajorVersion)

install(FILES
	${CMAKE_CURRENT_BINARY_DIR}/pumasConfig.cmake
	${CMAKE_CURRENT_BINARY_DIR}/pumasConfigVersion.cmake
	DESTINATION ${PUMAS_LIB}/cmake/pumas)

export(EXPORT pumasTarget
	FILE "${CMAKE_CURRENT_BINARY_DIR}/pumasTarget.cmake")

# Build and install rules for the examples, if enabled
if (PUMAS_BUILD_EXAMPLES)
        macro (pumas_example __tag __flux)
                set (__name "example-${__tag}")
                add_executable (${__name} "examples/pumas/${__tag}.c" ${__flux})
                target_compile_definitions (${__name} PRIVATE ${PUMAS_API})
                target_link_libraries (${__name} pumas)
                install (TARGETS ${__name} DESTINATION ${PUMAS_BIN} OPTIONAL)
        endmacro ()

        pumas_example (dump "")
        pumas_example (geometry "examples/pumas/flux.c")
        pumas_example (loader "")
        pumas_example (straight "examples/pumas/flux.c")
        pumas_example (tabulate "")
endif ()


# Build the test(s), if enabled
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        option (PUMAS_BUILD_TEST "Build the test suite" OFF)
        if (PUMAS_BUILD_TEST)
                include ("tests/test-pumas.cmake")
        endif ()
endif ()
