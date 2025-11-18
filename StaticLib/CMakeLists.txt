cmake_minimum_required(VERSION 3.14)

project(PlasmaLibrary LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (WIN32)
    set(lib_name Plasma_Windows_Lib)
endif (WIN32)
if (APPLE)
    set(lib_name Plasma_macOS_Lib)
endif (APPLE)
if (LINUX)
    set(lib_name Plasma_Linux_Lib)
endif (LINUX)

add_library(${lib_name} STATIC Plasma.hpp src/Plasma.cpp src/FFT.h src/FFT.cpp)
