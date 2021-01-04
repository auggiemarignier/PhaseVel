cmake_minimum_required(VERSION 3.10)

project(kernels)

set(CMAKE_VERBOSE_MAKEFILE ON)

enable_language(Fortran)

message("Fortran compiler is " ${CMAKE_Fortran_COMPILER})

set(MINEOS_DIR "" CACHE PATH "Path to MINEOS installation")
set(MINEOS_FDB_OBJ ${MINEOS_DIR}/fdb_eigen.o)
set(MINEOS_TIME_OBJ ${MINEOS_DIR}/time.o)

add_library(mineos_fdb OBJECT IMPORTED GLOBAL)
set_property(TARGET mineos_fdb PROPERTY IMPORTED_OBJECTS ${MINEOS_FDB_OBJ})
target_include_directories(mineos_fdb INTERFACE ${MINEOS_DIR}/fdb)

add_library(mineos_time OBJECT IMPORTED GLOBAL)
set_property(TARGET mineos_time PROPERTY IMPORTED_OBJECTS ${MINEOS_TIME_OBJ})
target_include_directories(mineos_time INTERFACE ${MINEOS_DIR}/time)

add_executable(${PROJECT_NAME} 
                src/kernels.f
                src/minos_bran.f
                src/eigcon_subs.f
                $<TARGET_OBJECTS:mineos_fdb>
                $<TARGET_OBJECTS:mineos_time>)

target_include_directories(${PROJECT_NAME} PUBLIC
    $<TARGET_PROPERTY:mineos_fdb,INTERFACE_INCLUDE_DIRECTORIES>
    $<TARGET_PROPERTY:mineos_time,INTERFACE_INCLUDE_DIRECTORIES>)
