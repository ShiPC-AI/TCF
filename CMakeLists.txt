# CMake
cmake_minimum_required(VERSION 3.0)
project(tcf)

# C++
set(CMAKE_CXX_STANDARD 17)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(-std=c++17 -g)
set(CMAKE_BUILD_TYPE "Release")

# multiple PCL exists, we select pcl-1.12
find_package(PCL 1.12 REQUIRED)
set(PCL_INCLUDE_DIRS /usr/local/include/pcl-1.12)
# find_package(PCL REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

# nlohmann_json
find_package(nlohmann_json REQUIRED)

add_subdirectory(utils)
add_subdirectory(registration)

add_executable(demo demo.cpp)
target_link_libraries(demo PRIVATE ${PCL_LIBRARIES} 
    ${TCF_REGISTRATION_LIB} ${TCF_UTILS_LIB} nlohmann_json::nlohmann_json) 