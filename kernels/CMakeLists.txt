cmake_minimum_required(VERSION 3.10)

project(kernels Fortran)
cmake_policy(SET CMP0097 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_Fortran_FLAGS -fallow-argument-mismatch)

include(ExternalProject)
ExternalProject_Add(mineos
    GIT_REPOSITORY git@github.com:geodynamics/mineos.git
    GIT_SHALLOW ON
    GIT_SUBMODULES ""

    CONFIGURE_COMMAND autoreconf -i <SOURCE_DIR> && <SOURCE_DIR>/configure FFLAGS=${CMAKE_Fortran_FLAGS} --disable-doc 
)