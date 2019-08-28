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
    EXPECT_EQ(expected, p3);

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
    std::string expected = "A[CCA,CCT,CG,T]T";
    EXPECT_EQ(expected,p);
}

/**
 * Two allelic classes segregate, each containing two SNPs
 * The graph with nodes of sequence length 1 has cross-class 'node sharing'
 * GOAL: collapse the node sharing with nested output showing the segregating classes.
 */
TEST(MSA, twoSegregatingClasses) {
    std::string MSA_string = ">R1\n"
                             "AAAAAAAAA\n"
                             ">R2\n"
                             "AATAAAAAA\n"
                             ">R3\n"
                             "AAAAATAAA\n"
                             ">R4\n"
                             "TTTTTTTTT\n"
                             ">R5\n"
                             "TTATTTTTT\n"
                             ">R6\n"
                             "TTTTTATTT\n";
    auto p = make_and_print_prg_string(MSA_string, false);
}

/**
 * An indel bubble ends outside of what, in its absence, would be another site.
 */
TEST(MSA, leakyBubble) {
    std::string MSA_string = ">R1\n"
                             "ATGCAT\n"
                             ">R2\n"
                             "AGCTAT\n"
                             ">R3\n"
                             "AG---T\n";
    auto p = make_and_print_prg_string(MSA_string, false);
}

/**
 * A set of indels all converging to one point
 */
TEST(MSA, expanded_IndelGraph) {
    std::string MSA_string = ">R1\n"
                             "AAAAAAAG\n"
                             ">R2\n"
                             "A------G\n"
                             ">R3\n"
                             "AA-----G\n"
                             ">R4\n"
                             "AAA----G\n"
                             ">R5\n"
                             "AAAA---G\n"
                             ">R6\n"
                             "AAAAA--G\n"
                             ">R7\n"
                             "AAAAAA-G\n";
    auto p = make_and_print_prg_string(MSA_string, false);
}

/**
 * For a low number of max incident bubbles parameter, the T/C SNP can get written several times.
 * GOAL: avoid T/C SNP getting written out several times
 */
TEST(MSA, SNP_rewriting) {
    std::string MSA_string = ">R1\n"
                             "GCCCA\n"
                             ">R2\n"
                             "G-CTA\n"
                             ">R3\n"
                             "G--CA\n";
    auto p = make_and_print_prg_string(MSA_string, false);
}

/**
 * A case where empty deletion record could occur.
 * GOAL: avoid empty allele
 */
TEST(MSA, successiveDeletions) {
    std::string MSA_string = ">R1\n"
                             "CGAACAAAG\n"
                             ">R2\n"
                             "CG--C--AG\n";
    auto p = make_and_print_prg_string(MSA_string, false);
 }

 /**
 * The reason we are doing all this in the first place.
 */
TEST(MSA, deletion_spanning_SNPs) {
    std::string MSA_string = ">R1\n"
                             "CGAACAAAG\n"
                             ">R2\n"
                             "CGTACATAG\n"
                             ">R3\n"
                             "CGAGCGAAG\n"
                             ">R4\n"
                             "CGAACCAAG\n"
                             ">R5\n"
                             "C-------G\n";
    auto p = make_and_print_prg_string(MSA_string, false);
}
