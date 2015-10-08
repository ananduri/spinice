cmake_minimum_required(VERSION 2.8)
project(kmc)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#Find Boost:
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIR})

#Find OpenMP:
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#Executables:
add_executable(a main.cpp flipsingle.cpp evolve.cpp evolvesave.cpp
getenergy.cpp mcstep.cpp mersenne.cpp)
target_link_libraries(a ${MKL_LIBRARIES})
