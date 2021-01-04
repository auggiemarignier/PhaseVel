cmake_minimum_required(VERSION 3.10)

project(kernels)

enable_language(Fortran)

message("Fortran compiler is " ${CMAKE_Fortran_COMPILER})

set(MINEOS_DIR "" CACHE PATH "Path to MINEOS installation")
set(MINEOS_OBJ ${MINEOS_DIR}/minos_bran.o)

add_library(mineos OBJECT IMPORTED GLOBAL)
set_property(TARGET mineos PROPERTY IMPORTED_OBJECTS ${MINEOS_OBJ})

add_executable(${PROJECT_NAME} src/kernels.f src/minos_bran.f)
