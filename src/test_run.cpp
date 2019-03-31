#include "msa_to_dfa.hpp"
#include <iostream>

int main(){

    MSA msa("/home/brice/Desktop/phD_Work/git_repos/prg_msa/test_data/AMA_Plasmodium.fasta");

    FA fa = FA(msa);


    for (auto node : fa.root.next){
        std::cout << node->letter << std::endl;
    }
    auto prg = oneDepth_prg(fa.root);
}