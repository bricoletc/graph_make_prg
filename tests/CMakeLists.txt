set(INCLUDE
        ${PROJECT_SOURCE_DIR}/inc)
add_executable(test_main
        prgs.cpp)
target_include_directories(test_main PUBLIC ${INCLUDE})
target_link_libraries(test_main libprg_msa)