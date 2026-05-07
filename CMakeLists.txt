cmake_minimum_required(VERSION 4.2)
project(untitled)

set(CMAKE_CXX_STANDARD 20)

add_library(untitled STATIC
        probable.cpp
        gallium_oxide.cpp
        probable.h
        vec3.h)