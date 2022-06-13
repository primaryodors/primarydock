# I like 3.11+ of CMake, since it features fetchcontent
# function, which can include 3rd party dependency libraries
# from github if necessary (e.g. approval tests)
cmake_minimum_required (VERSION 3.14)  # PI :)

# Name of project, programming language,
project(molecule_test LANGUAGES CXX)

# What C++ standard to use; see
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Make sure not to use any non-portable C++ extensions of the compiler
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "--coverage")

# Include class library
add_subdirectory(classes)


# This is how to declare an application (must include a main function)
add_executable(molecule_test molecule_test.cpp)

# This is how to link a library to an application
target_link_libraries(molecule_test class_library)

set(CMAKE_BUILD_TYPE RelWithDebInfo)

