include(${PROJECT_SOURCE_DIR}/ext/gtest.cmake)

set(INCLUDE
        ${PROJECT_SOURCE_DIR}/inc
        ./*.hpp
        ${source_dir}/googletest/include
        ${source_dir}/googlemock/include)

# To see gtest headers in IDE
include_directories(${INCLUDE})

file(GLOB_RECURSE SOURCE_FILES *.cpp)
add_executable(test_main
        ${SOURCE_FILES})
add_dependencies(test_main
        gtest)
target_include_directories(test_main PUBLIC
        ${INCLUDE})

target_link_libraries(test_main
        libprg_msa
        libgmock
        libgtest)