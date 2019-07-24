set (CMAKE_C_FLAGS_TEST
        "${CMAKE_C_FLAGS_DEBUG} --coverage" CACHE STRING
        "Flags used by the C compiler during test builds." FORCE)
mark_as_advanced(CMAKE_C_FLAGS_TEST)

include ("tests/libcheck.cmake")
include ("tests/materials.cmake")

add_executable (test-pumas EXCLUDE_FROM_ALL "tests/test-pumas.c")
target_compile_definitions (test-pumas PRIVATE -DPUMAS_VERSION=${PUMAS_VERSION})
add_dependencies (test-pumas LibCheck)
target_link_libraries (test-pumas check pumas)

if (${CMAKE_BUILD_TYPE} MATCHES Test)
        set (__args "coverage")
else ()
        set (__args "")
endif ()

cmake_policy (SET CMP0037 OLD)
if (${CMAKE_BUILD_TYPE} MATCHES "Test")
    # Add a target for test(s)
    add_custom_target (test DEPENDS test-pumas
        COMMAND lcov --directory . --zerocounters
        COMMAND test-pumas
        COMMAND lcov --directory . --capture --output-file coverage.info
        COMMAND lcov --remove coverage.info '*/tests/*' '*/include/*' '/usr/*' --output-file coverage.info
        COMMAND lcov --list coverage.info

        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running test(s)")
else ()
    add_custom_target (test DEPENDS test-pumas
        COMMAND test-pumas ${__args} 

        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running test(s)")
endif ()

add_custom_target (memcheck DEPENDS test-pumas
        COMMAND CK_NOFORK=on valgrind --error-exitcode=1 --leak-check=full --show-reachable=yes ./test-pumas

    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running memory check")
