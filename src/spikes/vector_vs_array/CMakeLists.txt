cmake_minimum_required (VERSION 3.14)
project(vector_vs_array LANGUAGES CXX)

# This is the main
add_library(library functions.cpp functions.hpp)
add_executable(vector_vs_array main.cpp)
target_link_libraries(vector_vs_array library)