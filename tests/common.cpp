#include "common.hpp"

std::vector<uint64_t> make_and_print_prg_string(const std::string &MSA_string, bool is_file) {
    MSA msa(MSA_string, is_file);
    FA fa(msa);
    sequence_Graph p = sequence_Graph(fa.root, MSA_string, is_file, 2);
    coverage_Graph p2 = coverage_Graph(p);
    // stringified_PRG<sequence_Graph, auto_Node> p3 = stringified_PRG<sequence_Graph, auto_Node>(p);
    stringified_PRG<coverage_Graph, coverage_Node> p4 = stringified_PRG<coverage_Graph, coverage_Node>(p2);
    auto& output = p4.prg;
    BOOST_LOG_TRIVIAL(info) << "Produced serialised PRG string of length " << output.size();
    if (output.size() < 100000){
        auto converted_to_string = to_readable_string(output);
        std::cout << converted_to_string << std::endl;
    }
    return output;
}
