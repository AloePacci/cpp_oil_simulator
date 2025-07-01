# Set the project name
project (map)
cmake_minimum_required(VERSION 3.10)

add_library(${PROJECT_NAME} map.cpp)
add_library(utils::map ALIAS ${PROJECT_NAME})

target_include_directories( ${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR} 
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Eigen
)
find_package(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Eigen)

target_link_libraries(${PROJECT_NAME}
    Eigen3::Eigen
)

file(GLOB MY_MAPS "${CMAKE_CURRENT_SOURCE_DIR}/../../mapas/*")

file(COPY ${MY_MAPS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../mapas)