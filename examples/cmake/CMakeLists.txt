# This example illustrates how to include PUMAS as a CMake package. Note that
# PUMAS must first have been installed using CMake, and that the install
# location must be known by CMake (e.g. by exporting it to CMAKE_PREFIX_PATH).
#
cmake_minimum_required (VERSION 2.8)
project (PUMAS_EXAMPLE C)

# Import the Pumas library (which must have been installed with CMake, first)
find_package(pumas CONFIG REQUIRED)

# Build the basic `tabulate` example
add_executable(example-tabulate ../pumas/tabulate.c)
target_link_libraries(example-tabulate pumas)
