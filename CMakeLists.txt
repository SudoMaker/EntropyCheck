cmake_minimum_required(VERSION 3.2)
project(EntropyCheck)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

add_library(EntropyCheck STATIC Source/EntropyCheck.c Source/EntropyCheck.h Source/EntropyCheck.hpp)
target_include_directories(EntropyCheck PUBLIC Source)

add_library(EntropyCheck_Shared SHARED Source/EntropyCheck.c Source/EntropyCheck.h Source/EntropyCheck.hpp)
set_target_properties(EntropyCheck_Shared PROPERTIES OUTPUT_NAME "EntropyCheck")

set(EntropyCheck_BuildTests -1)

if(${EntropyCheck_BuildTests} EQUAL -1)
    message("-- EntropyCheck: Building tests by default. Set EntropyCheck_BuildTests to 0 to disable.")
    set(EntropyCheck_BuildTests 1)
endif()

if(${EntropyCheck_BuildTests} EQUAL 1)
    add_executable(EntropyCheck_TestCpp Test/test_cpp.cpp)
    target_link_libraries(EntropyCheck_TestCpp EntropyCheck)
endif()