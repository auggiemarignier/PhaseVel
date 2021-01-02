cmake_minimum_required(VERSION 3.10)

project(Kernels)

enable_language(Fortran)

message(${CMAKE_Fortran_COMPILER})

set(MINEOS_DIR "" CACHE PATH "Path to MINEOS installation")


message(${MINEOS_DIR})