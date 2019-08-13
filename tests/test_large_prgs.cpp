#include "common.hpp"
#include "gtest/gtest.h"

TEST(MSA, Plasmodium_TwoRecords){
    //auto data_path = filesystem::path(__FILE__);
    std::string MSA_fpath = "/home/brice/Desktop/git_repos/prg_string_construction/prg_msa/tests/test_data/debugs/AMA_Plasmodium_cycle.fasta";
    make_and_print_prg_string(MSA_fpath);
}
