cmake_minimum_required(VERSION 3.10)

project(kernels Fortran)
cmake_policy(SET CMP0097 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_Fortran_FLAGS -fallow-argument-mismatch)

set(EP_MINEOS ep_mineos)
set(MINEOS mineos)

include(ExternalProject)
ExternalProject_Add(${EP_MINEOS}
    GIT_REPOSITORY https://github.com/geodynamics/mineos.git
    GIT_SHALLOW ON
    GIT_SUBMODULES ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND autoreconf -i <SOURCE_DIR> && <SOURCE_DIR>/configure FFLAGS=${CMAKE_Fortran_FLAGS} --disable-doc 
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(${EP_MINEOS} SOURCE_DIR)
include_directories(${SOURCE_DIR}/fdb)
ExternalProject_Get_Property(${EP_MINEOS} BINARY_DIR)
add_library(${MINEOS} OBJECT IMPORTED GLOBAL)
set_property(TARGET ${MINEOS} PROPERTY IMPORTED_OBJECTS ${BINARY_DIR}/fdb_eigen.o ${BINARY_DIR}/time.o)
add_dependencies(${MINEOS} ${EP_MINEOS})

add_executable(${PROJECT_NAME} 
                src/all_kernels.f
                src/kernelfunctions.f
                $<TARGET_OBJECTS:${MINEOS}>)
add_dependencies(${PROJECT_NAME} ${MINEOS})
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

install(TARGETS ${PROJECT_NAME})