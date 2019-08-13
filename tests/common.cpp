#include "common.hpp"
#include "oneDepth_prg.hpp"

std::string make_and_print_prg_string(const std::string &MSA_string, bool is_file) {
    MSA msa(MSA_string, is_file);
    FA fa(msa);
    nested_prg p = nested_prg(fa.root);
    auto& output = p.prg;
    BOOST_LOG_TRIVIAL(info) << "Produced PRG string of length " << output.length();
    if (output.length() < 10000){
        std::cout << p.prg << std::endl;
        std::cout << p.serialised_prg << std::endl;
    }
    return p.prg;
}
