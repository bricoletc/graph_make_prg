include(${PROJECT_SOURCE_DIR}/ext/gtest.cmake)

set(INCLUDE
        ${PROJECT_SOURCE_DIR}/inc
        ./*.hpp
        ${source_dir}/googletest/include
        ${source_dir}/googlemock/include)

# To see gtest headers in IDE
include_directories(${INCLUDE})

file(GLOB_RECURSE SOURCE_FILES ./*.cpp)
add_executable(test_make_prg
        ${SOURCE_FILES})
add_dependencies(test_make_prg
        gtest)
target_include_directories(test_make_prg PUBLIC
        ${INCLUDE})

target_link_libraries(test_make_prg
        lmake_prg
        libgmock
        libgtest)