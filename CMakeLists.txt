include(ExternalProject)
cmake_minimum_required(VERSION 3.13)
project(prg_msa)

# Right now I'm taking the approach of requiring Boost non head-only libraries to exist system-wide
# Down the line, downloading and compiling those locally to the project, as gramtools does, is preferable.

find_package(Threads REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.62.0 COMPONENTS filesystem log thread system REQUIRED)

file(GLOB_RECURSE SOURCE_FILES
        src/*.cpp)
set(INCLUDE ${PROJECT_SOURCE_DIR}/inc)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_SYSTEM_NO_DEPRECATED -fopenmp -Wall -Wextra")

# Make executable program
add_executable(prg_msa
    ${SOURCE_FILES})
target_link_libraries(prg_msa
        ${CMAKE_THREAD_LIBS_INIT}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_LOG_LIBRARY}
        ${Boost_THREAD_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${CMAKE_DL_LIBS}
        )
target_include_directories(prg_msa PUBLIC ${INCLUDE})

# Make linkable library
add_library(libprg_msa STATIC ${SOURCE_FILES})
target_include_directories(libprg_msa PUBLIC ${INCLUDE})
target_link_libraries(libprg_msa
        ${CMAKE_THREAD_LIBS_INIT}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_LOG_LIBRARY}
        ${Boost_THREAD_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${CMAKE_DL_LIBS}
        )

#Tests
enable_testing()
add_subdirectory(tests)
#add_test(test_main test_main)
