cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(FLFE VERSION 2.1.2 LANGUAGES C CXX)

if (PROJECT_BINARY_DIR STREQUAL PROJECT_SOURCE_DIR)
    message(WARNING "The binary directory of CMake cannot be the same as source directory!")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if (WIN32)
    add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES)
endif()

# without eigen
add_executable(main-mine)
file(GLOB_RECURSE sources CONFIGURE_DEPENDS src-mine/*.cpp src-mine/*.h)
target_sources(main-mine PUBLIC ${sources})

# with eigen
include_directories(~/miniconda3/envs/hdfr-stable/include/eigen3)
file(GLOB_RECURSE sources CONFIGURE_DEPENDS src-eigen/*.cpp src-eigen/*.h)
add_executable(main-eigen ${sources})
