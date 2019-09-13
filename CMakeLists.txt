include(ExternalProject)
cmake_minimum_required(VERSION 3.13)
project(graph_make_prg)

execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/download)
execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/src)

# Depend on sdsl
include(${PROJECT_SOURCE_DIR}/ext/sdsl.cmake)

# Right now I'm taking the approach of requiring Boost non head-only libraries to exist system-wide
# Down the line, downloading and compiling those locally to the project, as gramtools does, is preferable.

find_package(Threads REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.62.0 COMPONENTS filesystem log thread system REQUIRED)

file(GLOB_RECURSE SOURCE_FILES
        src/*.cpp)
set(INCLUDE ${PROJECT_SOURCE_DIR}/inc)
set(EXTERNAL_INCLUDE_DIR
        ${CMAKE_CURRENT_BINARY_DIR}/include)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_SYSTEM_NO_DEPRECATED -fopenmp -Wall -Wextra")

# Make linkable library
add_library(lmake_prg STATIC ${SOURCE_FILES})
target_include_directories(lmake_prg PUBLIC ${INCLUDE} ${EXTERNAL_INCLUDE_DIR})
add_dependencies(lmake_prg
        sdsl)
target_link_libraries(lmake_prg
        ${CMAKE_THREAD_LIBS_INIT}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_LOG_LIBRARY}
        ${Boost_THREAD_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${CMAKE_DL_LIBS}
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libsdsl.a
        )

# Make executable program
add_executable(make_prg
        ${SOURCE_FILES})
target_link_libraries(make_prg
        ${CMAKE_THREAD_LIBS_INIT}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_LOG_LIBRARY}
        ${Boost_THREAD_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${CMAKE_DL_LIBS}
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libsdsl.a
        )
target_include_directories(make_prg PUBLIC ${INCLUDE})

#Tests
enable_testing()
add_subdirectory(tests)
#add_test(test_main test_main)
