include(ExternalProject)

set (CHECK_DIR "${CMAKE_BINARY_DIR}/check")

ExternalProject_Add (LibCheck
    GIT_REPOSITORY "https://github.com/libcheck/check"
    GIT_TAG "0.12.0"
    INSTALL_DIR ${CHECK_DIR}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CHECK_DIR}
)

# Prevent a bug in CMake raising an error because the include dir does not yet
# exist (see: https://cmake.org/Bug/view.php?id=15052)
file (MAKE_DIRECTORY "${CHECK_DIR}/include")

add_library (check STATIC IMPORTED)
set_target_properties (check PROPERTIES
    IMPORTED_LOCATION "${CHECK_DIR}/lib/libcheck.a"
    INTERFACE_INCLUDE_DIRECTORIES "${CHECK_DIR}/include"
    INTERFACE_LINK_LIBRARIES "-lm -lrt"
)
