include(ExternalProject)
cmake_minimum_required(VERSION 3.13)
project(prg_msa)

set(CMAKE_CXX_STANDARD 14)

file(GLOB_RECURSE SOURCE_FILES
        src/*.cpp)

set(INCLUDE inc)

add_executable(prg_msa
    ${SOURCE_FILES})
target_include_directories(prg_msa PUBLIC ${INCLUDE})

add_library(libprg_msa ${SOURCE_FILES})
target_include_directories(libprg_msa PUBLIC ${INCLUDE})

#Tests
enable_testing()
add_subdirectory(tests)
#add_test(test_main test_main)
