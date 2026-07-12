include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

add_executable(geometry_kernel_unit_tests
    core/test_types.cpp
)
target_link_libraries(geometry_kernel_unit_tests PRIVATE GTest::gtest_main geometry_kernel)
target_include_directories(geometry_kernel_unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
gtest_discover_tests(geometry_kernel_unit_tests)
