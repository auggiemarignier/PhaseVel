cmake_minimum_required(VERSION 3.10)

project(kernels)
cmake_policy(SET CMP0097 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)

enable_language(Fortran)

include(ExternalProject)
ExternalProject_Add(mineos
    GIT_REPOSITORY git@github.com:geodynamics/mineos.git
    GIT_SHALLOW ON
    GIT_SUBMODULES ""
)