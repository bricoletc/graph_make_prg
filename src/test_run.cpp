#include "msa_to_dfa.hpp"
#include <iostream>

int main(){

    MSA msa("/home/brice/Desktop/git_repos/prg_msa/test_data/AMA_Plasmodium.fasta");

    auto root = MSA_to_FA(msa);
}