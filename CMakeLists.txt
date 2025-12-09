cmake_minimum_required(VERSION 2.8)

project(pathPlanningAlgorithm)

#json库
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../common/json
    ${CMAKE_CURRENT_SOURCE_DIR}/../../common/logger/inc
    ${CMAKE_CURRENT_SOURCE_DIR}/../../common/logger
    ${CMAKE_CURRENT_SOURCE_DIR}/../../common/config)

find_package(Eigen3 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(GDAL)

add_library(algo uavPathPlanning.cpp)

add_library(math_util 
    math_util/polygon2d.cpp
    math_util/box2d.cpp
    math_util/linear_interpolation.cpp
    math_util/math_utils.cpp
    math_util/minimum_snap.cpp
    )

# minimum_snap.cpp now contains trajectory-generator and planner implementations
target_sources(math_util PRIVATE math_util/minimum_snap.cpp)

target_link_libraries(algo math_util)
target_link_libraries(math_util yaml-cpp)
if (GDAL_FOUND)
    target_include_directories(algo PUBLIC ${GDAL_INCLUDE_DIRS})
    target_link_libraries(algo ${GDAL_LIBRARIES})
    add_definitions(-DHAVE_GDAL)
endif()

if (GDAL_FOUND)
    message(STATUS "GDAL found, but no dedicated elevation_to_costmap tool source present in this tree. If you add one at math_util/elevation_to_costmap.cpp it will be built.")
else()
    message(STATUS "GDAL not found: GDAL-based tools will be unavailable")
endif()

add_executable(uavPathPlanningTest uavPathPlanningTest.cpp)

# target_link_libraries(uavPathPlanningTest 
#     algo
#     Eigen3::Eigen
#     nlohmann_json::nlohmann_json
# )
target_link_libraries(uavPathPlanningTest algo)
