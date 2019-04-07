#include "MSA.hpp"
#include "FA.hpp"
#include "oneDepth_prg.hpp"
#include "nested_prg.hpp"
#include <iostream>


int main(){

    //MSA msa("/home/brice/Desktop/phD_Work/git_repos/prg_msa/test_data/AMA_Plasmodium.fasta");
    //MSA msa("/home/brice/Desktop/phD_Work/git_repos/prg_msa/test_data/AMA_Plasmodium_2recs.fasta");
    //MSA msa("/home/brice/Desktop/phD_Work/git_repos/prg_msa/test_data/nested.fasta");

    //MSA msa("/home/brice/Desktop/git_repos/prg_msa/test_data/bifurc.fasta");
    //MSA msa("/home/brice/Desktop/git_repos/prg_msa/test_data/all_ends.fasta");
    //MSA msa("/home/brice/Desktop/phD_Work/git_repos/prg_msa/test_data/bifurc2.fasta");
    MSA msa("/home/brice/Desktop/phD_Work/git_repos/prg_msa/test_data/bifurc3.fasta");

    FA fa = FA(msa);

    nested_prg p = nested_prg(fa.root);
    //auto prg = oneDepth_prg(fa.root);
    //std::cout << prg.prg << std::endl;

    std::cout << "Final prg: " << p.prg << std::endl;
    //void deletion_above_snpAndIndel();
    //deletion_above_snpAndIndel();

}

#include <sstream>

void deletion_above_snpAndIndel(){
    std::string MSA_string = ">Rec1\n"
                             "ACGTTA\n"
                             ">Rec2\n"
                             "ACC-TA\n"
                             ">Rec3\n"
                             "A----A\n"
    ;
    std::string test_DIR = "/home/brice/Desktop/git_repos/prg_msa/test_data/";


    std::string tmp_file = test_DIR + "tmp.fasta";

    std::ofstream tmp_out = std::ofstream(tmp_file);
    tmp_out << MSA_string;
    tmp_out.close();

    MSA msa(tmp_file);

    FA fa = FA(msa);

    auto prg = oneDepth_prg(fa.root);
    std::cout << prg.prg << std::endl;
}
