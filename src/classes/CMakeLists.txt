# I like 3.11+ of CMake, since it features fetchcontent
# function, which can include 3rd party dependency libraries
# from github if needed (e.g. approval tests)
cmake_minimum_required (VERSION 3.14)  # PI :)

# Name of project, programming language,
project(class_library LANGUAGES CXX)

# What C++ standard to use; see
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Make sure not to use any non-portable C++ extensions of the compiler
set(CMAKE_CXX_EXTENSIONS OFF)

# This is to enable code coverage analysis in CLion
set(CMAKE_CXX_FLAGS "--coverage")

# Assign all relevant source files to variable
file(GLOB library_source *.h *.cpp)

# This is how to declare a library
add_library(class_library ${library_source})

set(CMAKE_BUILD_TYPE RelWithDebInfo)

