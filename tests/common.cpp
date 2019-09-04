#include "common.hpp"

std::vector<uint64_t> make_and_print_prg_string(const std::string &MSA_string, bool is_file) {
    MSA msa(MSA_string, is_file);
    FA fa(msa);
    sequence_Graph p = sequence_Graph(fa.root, MSA_string, is_file, 2);
    coverage_Graph p2 = coverage_Graph(p);
    stringified_PRG<sequence_Graph> p3 = stringified_PRG<sequence_Graph>(p);
    auto& output = p3.prg;
    BOOST_LOG_TRIVIAL(info) << "Produced unserialised PRG string of length " << output.size();
    if (output.size() < 100000){
        std::cout << output.size() << std::endl;
    }
    return output;
}
