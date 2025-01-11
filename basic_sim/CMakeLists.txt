cmake_minimum_required(VERSION 3.10)
project(
  BasicSimulator
  VERSION 1.0
  LANGUAGES CXX)

# target_include_directories(Tutorial PUBLIC
#                            "${PROJECT_BINARY_DIR}"
#                            )
add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME}
    utils::map
    utils::basic_simulator
)
# target_include_directories(basic_sim ../${CMAKE_CURRENT_SOURCE_DIR})
