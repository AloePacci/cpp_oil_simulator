# Set the project name
project (basic_simulator)


add_library(${PROJECT_NAME} basic_simulator.cpp)
add_library(utils::basic_simulator ALIAS ${PROJECT_NAME})
find_package(Eigen3 REQUIRED)

target_include_directories( ${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}
)

target_link_libraries(${PROJECT_NAME}
    utils::map
    Eigen3::Eigen
)