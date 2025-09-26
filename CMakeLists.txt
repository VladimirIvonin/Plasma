cmake_minimum_required(VERSION 3.14)

project(Plasma LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(Plasma SHARED Plasma_funcs.hpp Plasma_funcs.cpp FFT.h FFT.cpp)

target_compile_definitions(Plasma PRIVATE PLASMA_LIBRARY)


if (WIN32)
    set(CMAKE_CXX_FLAGS "-static-libgcc -static-libstdc++")
endif()
