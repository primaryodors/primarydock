cmake_minimum_required (VERSION 3.14)

# Name of project, programming language,
project(vector_vs_array LANGUAGES CXX)

# What C++ standard to use; see
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Make sure not to use any non-portable C++ extensions of the compiler
set(CMAKE_CXX_EXTENSIONS OFF)

# This is how to declare a library
add_library(library functions.cpp functions.hpp)

# This is how to declare an application (must include a main function)
add_executable(vector_vs_array main.cpp)

# This is how to link a library to an application
target_link_libraries(vector_vs_array library)