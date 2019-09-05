#include "common.hpp"
#include "gtest/gtest.h"

std::pair<std::string, std::string> make_and_print_prg_string(const std::string &MSA_string, bool is_file) {
    MSA msa(MSA_string, is_file);
    FA fa(msa);
    sequence_Graph p = sequence_Graph(fa.root, MSA_string, is_file, 2);
    coverage_Graph p2 = coverage_Graph(p);
    stringified_PRG<sequence_Graph, auto_Node> p3 = stringified_PRG<sequence_Graph, auto_Node>(p);
    stringified_PRG<coverage_Graph, coverage_Node> p4 = stringified_PRG<coverage_Graph, coverage_Node>(p2);
    auto output_seqG = to_readable_string(p3.prg);
    auto output_covG = to_readable_string(p4.prg);
    BOOST_LOG_TRIVIAL(debug) << "Produced serialised PRG string of length " << output_covG.size();
    if (output_covG.size() < 100000){
        std::cout << output_covG << std::endl;
    }
    return std::make_pair(output_seqG, output_covG);
}
