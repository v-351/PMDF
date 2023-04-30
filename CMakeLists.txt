cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "mpicxx")
set( MPIEXEC_EXECUTABLE "mpirun" )

project(PMDF LANGUAGES CXX)
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED on )

include_directories(include)

find_package(MPI)

add_library(CELL SHARED src/cell.cpp)
add_library(LIB SHARED src/particle.cpp src/models.cpp src/property.cpp src/metrics.cpp src/mpi-wrapper.cpp)
add_library(MODELS SHARED src/models.cpp)
add_library(STL_SURFACE SHARED src/surface.cpp)
add_executable( ${PROJECT_NAME} src/main.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC CELL)
target_link_libraries(${PROJECT_NAME} PUBLIC LIB)
target_link_libraries(${PROJECT_NAME} PUBLIC MODELS)
target_link_libraries(${PROJECT_NAME} PUBLIC STL_SURFACE)
target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)