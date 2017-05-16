
cmake_minimum_required (VERSION 2.6)
project (Pumas)

# Build options.
set (BUILD_SHARED_LIBS TRUE
    CACHE BOOL "Build the library as shared")
option (BUILD_TABULATE
    "Build the pumas-tabulate utility" ON)
if (BUILD_TABULATE)
    add_definitions(-D_BUILD_TABULATE)
endif (BUILD_TABULATE)
if ("${CMAKE_SYSTEM}" MATCHES "Linux")
    option(GDB_MODE
    "Additional features for debugging with gdb" OFF)
    if (GDB_MODE)
        add_definitions(-DGDB_MODE)
    endif (GDB_MODE)
endif ("${CMAKE_SYSTEM}" MATCHES "Linux")

# Build and install rules for the PUMAS library.
include_directories (include)
add_library (pumas src/pumas.c include/pumas.h)
if (WIN32)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    if (BUILD_SHARED_LIBS)
        target_compile_definitions(pumas
            PUBLIC "-DPUMAS_API=__declspec(dllexport)")
    endif (BUILD_SHARED_LIBS)
    install (TARGETS pumas DESTINATION .)
    install (FILES include/pumas.h DESTINATION .)
else (WIN32)
    install (TARGETS pumas DESTINATION lib)
    install (FILES include/pumas.h DESTINATION include)
endif (WIN32)

# Build and install rules for the tabulate utility, if enabled.
if (BUILD_TABULATE)
    add_executable (pumas-tabulate src/pumas-tabulate.c src/pumas-tabulate.h)
    if (WIN32)
        if (BUILD_SHARED_LIBS)
            target_compile_definitions(pumas-tabulate
                PUBLIC "-DPUMAS_API=__declspec(dllimport)")
        endif (BUILD_SHARED_LIBS)
        target_link_libraries (pumas-tabulate pumas)
        install (TARGETS pumas-tabulate DESTINATION .)
    else (WIN32)
        target_link_libraries (pumas-tabulate m)
        target_link_libraries (pumas-tabulate pumas)
        install (TARGETS pumas-tabulate DESTINATION bin)
    endif (WIN32)
endif (BUILD_TABULATE)
