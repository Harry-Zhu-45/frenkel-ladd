cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(FLFE VERSION 3.0.0 LANGUAGES C CXX)

if (PROJECT_BINARY_DIR STREQUAL PROJECT_SOURCE_DIR)
    message(WARNING "The binary directory of CMake cannot be the same as source directory!")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if (WIN32)
    add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES)
endif()

# My Box
include_directories(~/.conda/envs/hoomd4-stable/include/eigen3)
file(GLOB_RECURSE sources CONFIGURE_DEPENDS src/MyBox.cpp)
add_executable(main-MyBox ${sources})

# Original Box
include_directories(~/.conda/envs/hoomd4-stable/include/eigen3)
file(GLOB_RECURSE sources CONFIGURE_DEPENDS src/OriginalBox.cpp)
add_executable(main-OriginalBox ${sources})
