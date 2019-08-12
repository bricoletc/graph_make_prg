#include <nested_prg.hpp>
#include "oneDepth_prg.hpp"

void deletion_above_snpAndIndel(){
    bool is_file = false;
    std::string MSA_string = ">Rec1\n"
                      "ACGTTA\n"
                      ">Rec2\n"
                      "ACC-TA\n"
                      ">Rec3\n"
                      "A----A\n"
                      ;

    MSA msa(MSA_string, is_file);
    FA fa(msa);
    nested_prg p = nested_prg(fa.root);
}

int main(){
    deletion_above_snpAndIndel();
}