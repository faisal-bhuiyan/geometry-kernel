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
    core/test_tolerance.cpp
    core/test_predicates2.cpp
    core/test_predicates3.cpp
    core/test_area.cpp
    core/test_point2.cpp
    core/test_point3.cpp
    core/test_vector2.cpp
    core/test_vector3.cpp
    queries/test_closest_point2.cpp
    queries/test_closest_point3.cpp
    queries/test_containment2.cpp
    queries/test_containment3.cpp
    queries/test_intersection.cpp
    algorithms/test_polygon_mesh.cpp
    algorithms/test_ear_clipping.cpp
)
target_link_libraries(geometry_kernel_unit_tests PRIVATE GTest::gtest_main geometry_kernel)
target_include_directories(geometry_kernel_unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
gtest_discover_tests(geometry_kernel_unit_tests)
