if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
        target_compile_options(pumas PUBLIC "--coverage")
        if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
                target_link_options(pumas PUBLIC --coverage)
        else()
                target_link_libraries(pumas PUBLIC --coverage)
        endif()
endif ()

include ("tests/libcheck.cmake")

add_executable (test-pumas EXCLUDE_FROM_ALL "tests/test-pumas.c")
target_compile_definitions (test-pumas PRIVATE
        -DPUMAS_VERSION_MAJOR=${PUMAS_VERSION_MAJOR}
        -DPUMAS_VERSION_MINOR=${PUMAS_VERSION_MINOR}
        -DPUMAS_VERSION_PATCH=${PUMAS_VERSION_PATCH})
add_dependencies (test-pumas LibCheck)
target_link_libraries (test-pumas check pumas)

add_custom_command(
        TARGET test-pumas POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/data/materials.xml
                ${CMAKE_CURRENT_BINARY_DIR}/materials/materials.xml)

add_custom_command(
        TARGET test-pumas POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
                ${CMAKE_CURRENT_BINARY_DIR}/materials/dedx/muon)

add_custom_command(
        TARGET test-pumas POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
                ${CMAKE_CURRENT_BINARY_DIR}/materials/dedx/tau)

if (${CMAKE_VERSION} VERSION_GREATER "3.0.0")
        cmake_policy (SET CMP0037 OLD)
endif ()

# Add a target for test(s)
add_custom_target (test DEPENDS test-pumas
        COMMAND test-pumas

        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running test(s)")

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
        # Add a target for coverage (requires lcov)
        add_custom_target (coverage DEPENDS test-pumas
                COMMAND lcov --directory . --zerocounters
                COMMAND test-pumas
                COMMAND lcov --directory . --capture --output-file coverage.info
                COMMAND lcov --remove coverage.info '*/tests/*' '*/include/*' '/usr/*' --output-file coverage.info
                COMMAND lcov --list coverage.info

                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                COMMENT "Running coverage")
endif ()

# Add a target for testing memory usage (requires valgrind)
add_custom_target (memcheck DEPENDS test-pumas
        COMMAND CK_NOFORK=on valgrind --error-exitcode=1 --leak-check=full --show-reachable=yes ./test-pumas

        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running memory check")
