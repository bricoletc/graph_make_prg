#include <nested_prg.hpp>
#include <assert.h>
#include "oneDepth_prg.hpp"
#include "log.hpp"

//TODO: an evaluator of PRG string correctness which is independent of allele ordering- or, order the alleles alphabetically.

std::string make_and_print_prg_string(const std::string& MSA_string){
    bool is_file = false;
    MSA msa(MSA_string, is_file);
    FA fa(msa);
    nested_prg p = nested_prg(fa.root);
    std::cout << p.prg << std::endl;
    return p.prg;
}


/*
 * Testing a series here, of related but slightly different NFAs
 */
void deletion_spanning_snpAndIndel(){
    std::string MSA_string = ">Rec1\n"
                      "ACGTTA\n"
                      ">Rec2\n"
                      "ACC-TA\n"
                      ">Rec3\n"
                      "A----A\n"
                      ;
    make_and_print_prg_string(MSA_string);

    // A single difference in Rec 3 which turns out to be important for ability to make correct PRG strings
    std::string MSA_string2 = ">Rec1\n"
                             "ACGTTA\n"
                             ">Rec2\n"
                             "ACC-TA\n"
                             ">Rec3\n"
                             "A---TA\n"
    ;
    make_and_print_prg_string(MSA_string2);


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

    auto p3 = make_and_print_prg_string(MSA_string3);
    std::string expected = "A[C[GT,C]T[A,G],A]";
    assert(expected == p3);

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
    init_logging("debug");
    deletion_spanning_snpAndIndel();
    adjacent_snp_and_del();
}