#include "common.hpp"
#include <assert.h>
#include "gtest/gtest.h"


/*
 * Testing a series here, of related but slightly different NFAs
 */
TEST(MSA, deletion_spanning_snpAndIndel) {
    std::string MSA_string = ">Rec1\n"
                             "ACGTTA\n"
                             ">Rec2\n"
                             "ACC-TA\n"
                             ">Rec3\n"
                             "A----A\n";
    make_and_print_prg_string(MSA_string, false);

    // A single difference in Rec 3 which turns out to be important for ability to make correct PRG strings
    std::string MSA_string2 = ">Rec1\n"
                              "ACGTTA\n"
                              ">Rec2\n"
                              "ACC-TA\n"
                              ">Rec3\n"
                              "A---TA\n";
    make_and_print_prg_string(MSA_string2, false);


    std::string MSA_string3 = ">Rec1\n"
                              "ACGTTA\n"
                              ">Rec2\n"
                              "ACC-TA\n"
                              ">Rec3\n"
                              "ACC-TG\n"
                              ">Rec4\n"
                              "A----A\n"
                              ">Rec5\n"
                              "ACGTTG\n";

    auto p3 = make_and_print_prg_string(MSA_string3, false);
    std::string expected = "A[A,C[C,GT]T[A,G]]";
    assert(expected == p3);

}

TEST(MSA, adjacent_snp_and_del) {
    std::string MSA_string = ">Rec1\n"
                             "ACGT\n"
                             ">Rec2\n"
                             "ACCT\n"
                             ">Rec3\n"
                             "A-CT\n";

    auto p3 = make_and_print_prg_string(MSA_string, false);
}

TEST(MSA, level3Nesting_allEndAtFinalNode) {

    std::string MSA_string = ">R1\n"
                             "ACCAT\n"
                             ">R2\n"
                             "ACG-T\n"
                             ">R3\n"
                             "AT--T\n"
                             ">R4\n"
                             "ACCTT\n";
    auto p = make_and_print_prg_string(MSA_string, false);
    std::string expected = "A[C[C[A,T],G],T]T";
    assert(expected == p);
}

/**
 * Two allelic classes segregate, each containing two SNPs
 * The graph with nodes of sequence length 1 looks horrible.
 */
TEST(MSA, twoSegregatingClasses) {
    std::string MSA_string = ">R1\n"
                             "ATCCGA\n"
                             ">R2\n"
                             "ATCCAA\n"
                             ">R3\n"
                             "ATCAGA\n"
                             ">R4\n"
                             "ACGGTA\n"
                             ">R5\n"
                             "ACGATA\n"
                             ">R6\n"
                             "ACGGAA\n";
    auto p = make_and_print_prg_string(MSA_string, false);
}
