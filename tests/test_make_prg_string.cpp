/**
 * @file Builds PRG graphs from MSAs: a sequence graph and a coverage graph from the sequence graph.
 * Stringifies both and tests for string equality.
 */
#include "common.hpp"
#include <assert.h>
#include "gtest/gtest.h"


/**
 * -----------------------
 * Small, in-string PRGs
 * -----------------------
 */

/*
 * Testing a series here, of related but slightly different NFAs
 */
TEST(MSA, deletion_spanning_snpAndIndel) {
    std::pair<std::string, std::string> res;
    std::string MSA_string = ">Rec1\n"
                             "ACGTTA\n"
                             ">Rec2\n"
                             "ACC-TA\n"
                             ">Rec3\n"
                             "A----A\n";
    res = make_and_print_prg_string(MSA_string, false);
    EXPECT_EQ(res.first, res.second);

    // A single difference in Rec 3 which turns out to be important for ability to make correct PRG strings
    std::string MSA_string2 = ">Rec1\n"
                              "ACGTTA\n"
                              ">Rec2\n"
                              "ACC-TA\n"
                              ">Rec3\n"
                              "A---TA\n";
    res = make_and_print_prg_string(MSA_string2, false);
    EXPECT_EQ(res.first, res.second);


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

    res = make_and_print_prg_string(MSA_string3, false);
    EXPECT_EQ(res.first, res.second);
    std::string expected = "A[A,C[C,GT]T[A,G]]";
    EXPECT_EQ(expected, res.second);

}

TEST(MSA, adjacent_snp_and_del) {
    std::string MSA_string = ">Rec1\n"
                             "ACGT\n"
                             ">Rec2\n"
                             "ACCT\n"
                             ">Rec3\n"
                             "A-CT\n";

    auto res = make_and_print_prg_string(MSA_string, false);
    EXPECT_EQ(res.first, res.second);
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
    auto res = make_and_print_prg_string(MSA_string, false);
    EXPECT_EQ(res.first, res.second);
    std::string expected = "A[CCA,CCT,CG,T]T";
    //EXPECT_EQ(expected,p);
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
    auto res = make_and_print_prg_string(MSA_string, false);
    EXPECT_EQ(res.first, res.second);
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
    auto res = make_and_print_prg_string(MSA_string, false);
    EXPECT_EQ(res.first, res.second);
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
    auto res = make_and_print_prg_string(MSA_string, false);
    EXPECT_EQ(res.first, res.second);
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
    auto res = make_and_print_prg_string(MSA_string, false);
    EXPECT_EQ(res.first, res.second);
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
    auto res = make_and_print_prg_string(MSA_string, false);
    EXPECT_EQ(res.first, res.second);
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
    auto res = make_and_print_prg_string(MSA_string, false);
    EXPECT_EQ(res.first, res.second);
}

/**
 * -----------------------
 * Large, in-file PRGs
 * -----------------------
 */

/**
 * The AMA1 gene of P falciparum
 */
TEST(MSA, Plasmodium_TwoRecords){
    //auto data_path = filesystem::path(__FILE__);
    //std::string MSA_fpath = "/home/brice/Desktop/git_repos/prg_string_construction/prg_msa/tests/test_data/AMA_Plasmodium_cycle.fasta";
    std::string MSA_fpath = "/home/brice/Desktop/git_repos/prg_string_construction/prg_msa/tests/test_data/AMA_Plasmodium.fasta";
    auto res = make_and_print_prg_string(MSA_fpath);
    EXPECT_EQ(res.first, res.second);
}
