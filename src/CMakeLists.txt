project(Surrender)
cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(PYBIND11_PYTHON_VERSION 3.6)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC")

add_subdirectory(utils)
add_subdirectory(basic_sim)