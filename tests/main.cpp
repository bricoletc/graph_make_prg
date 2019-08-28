#include "gtest/gtest.h"
#include "common.hpp"

int main(int argc, char **argv) {
    init_logging("info");
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
