cmake_minimum_required(VERSION 3.13)
project(prg_msa)

set(CMAKE_CXX_STANDARD 14)

include_directories(inc)

add_executable(prg_msa
    inc/msa_to_dfa.hpp
    src/msa_to_dfa.cpp
    src/test_run.cpp)
