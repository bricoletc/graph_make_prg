#include "common.hpp"

std::string make_and_print_prg_string(const std::string &MSA_string, bool is_file) {
    MSA msa(MSA_string, is_file);
    FA fa(msa);
    nested_prg p = nested_prg(fa.root);
    std::cout << p.prg << std::endl;
    return p.prg;
}
