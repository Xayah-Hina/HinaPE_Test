cmake_minimum_required(VERSION 3.26)
cmake_policy(SET CMP0091 NEW)
set(CMAKE_MSVC_RUNTIME_LIBRARY MultiThreadedDLL)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL /fp:fast /openmp /arch:AVX2") # MSVC
project(HinaPE_Test)
include(./FindHoudini.cmake)

# Find CUDA
find_package(CUDAToolkit REQUIRED)
enable_language(CUDA)
if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75)
endif ()
set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)
set(CUDA_SEPARABLE_COMPILATION ON)
set(CUDA_STANDARD 17)

add_subdirectory(CUDA_CubbyFlow)

add_library(
        SIM_TEST
        SHARED
        hdk/Entrance.cpp
        hdk/GAS_TestCubbyFlowSmoke.cpp
        hdk/GAS_TestCubbyFlowSmoke.h
        hdk/GAS_TestCubbyFlowFLIP.cpp
        hdk/GAS_TestCubbyFlowFLIP.h
        hdk/GAS_TestPhiFlow.cpp
        hdk/GAS_TestPhiFlow.h
        hdk/cubby.cpp
        hdk/cubby.h
)
target_link_libraries(
        SIM_TEST
        PUBLIC
        Houdini
        CUDA_CubbyFlow
)
target_include_directories(SIM_TEST PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

houdini_configure_target(SIM_TEST)
