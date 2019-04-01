#include "MSA.hpp"
#include "FA.hpp"
#include "oneDepth_prg.hpp"
#include <iostream>

int main(){

    //MSA msa("/home/brice/Desktop/git_repos/prg_msa/test_data/AMA_Plasmodium.fasta");
    //MSA msa("/home/brice/Desktop/git_repos/prg_msa/test_data/AMA_Plasmodium_2recs.fasta");
    //MSA msa("/home/brice/Desktop/git_repos/prg_msa/test_data/GC00001222_na_aln.fa");
     MSA msa("/home/brice/Desktop/git_repos/prg_msa/test_data/nested.fasta");


    FA fa = FA(msa);

    auto prg = oneDepth_prg(fa.root);
    std::cout << prg.prg << std::endl;
}