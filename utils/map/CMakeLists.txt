# Set the project name
project (map)


add_library(${PROJECT_NAME} map.cpp)
add_library(utils::map ALIAS ${PROJECT_NAME})

target_include_directories( ${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}
)
find_package(Eigen3 REQUIRED)

target_link_libraries(${PROJECT_NAME}
    Eigen3::Eigen
)

file(GLOB MY_MAPS "${CMAKE_CURRENT_SOURCE_DIR}/../../mapas/*")

file(COPY ${MY_MAPS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../mapas)