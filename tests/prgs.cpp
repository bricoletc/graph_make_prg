#include <nested_prg.hpp>
#include <assert.h>
#include "oneDepth_prg.hpp"

//TODO: an evaluator of PRG string correctness which is independent of allele ordering- or, order the alleles alphabetically.

/*
 * Testing a series here, of related but slightly different NFAs
 */
void deletion_spanning_snpAndIndel(){
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

    // A single difference in Rec 3 which turns out to be important for ability to make correct PRG strings
    std::string MSA_string2 = ">Rec1\n"
                             "ACGTTA\n"
                             ">Rec2\n"
                             "ACC-TA\n"
                             ">Rec3\n"
                             "A---TA\n"
    ;

    MSA msa2(MSA_string2, is_file);
    FA fa2(msa2);
    nested_prg p2 = nested_prg(fa2.root);

    std::string MSA_string3 = ">Rec1\n"
                              "ACGTTA\n"
                              ">Rec2\n"
                              "ACC-TA\n"
                              ">Rec3\n"
                              "ACC-TG\n"
                              ">Rec4\n"
                              "A----A\n"
                              ">Rec5\n"
                              "ACGTTG\n"
    ;
    MSA msa3(MSA_string3, is_file);
    FA fa3(msa3);
    nested_prg p3 = nested_prg(fa3.root);
    std::string expected = "A[C[GT,C]T[A,G],A]";
    assert(expected == p3.prg);

}

void adjacent_snp_and_del(){
    bool is_file = false;
    std::string MSA_string = ">Rec1\n"
                             "ACGT\n"
                             ">Rec2\n"
                             "ACCT\n"
                             ">Rec3\n"
                             "A-CT\n";

    MSA msa(MSA_string, is_file);
    FA fa(msa);
    nested_prg p = nested_prg(fa.root);
}

int main(){
    deletion_spanning_snpAndIndel();
    adjacent_snp_and_del();
}