# Set the project name
project (basic_simulator)
cmake_minimum_required(VERSION 3.10)
add_library(${PROJECT_NAME} basic_simulator.cpp)
add_library(utils::basic_simulator ALIAS ${PROJECT_NAME})
find_package(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Eigen)

target_include_directories( ${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Eigen
)

target_link_libraries(${PROJECT_NAME}
    utils::map
    Eigen3::Eigen
)