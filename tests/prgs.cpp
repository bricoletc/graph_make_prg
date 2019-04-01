#include "oneDepth_prg.hpp"

#include <sstream>

void deletion_above_snpAndIndel(){
    std::string MSA_string = ">Rec1\n"
                      "ACGTTA\n"
                      ">Rec2\n"
                      "ACC-TA\n"
                      ">Rec3\n"
                      "A----A\n"
                      ;

    std::stringstream msa_file = std::stringstream(MSA_string);
}