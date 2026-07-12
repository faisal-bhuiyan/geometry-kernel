cmake_minimum_required(VERSION 3.20)

project(GeometryKernel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_library(geometry_kernel INTERFACE)
target_include_directories(geometry_kernel INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src)

include(CTest)
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
