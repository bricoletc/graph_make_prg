cmake_minimum_required(VERSION 3.13)
project(prg_msa)

set(CMAKE_CXX_STANDARD 14)

include_directories(inc)

add_executable(prg_msa
        src/MSA.cpp
        src/FA.cpp
        src/oneDepth_prg.cpp
    src/test_run.cpp)
