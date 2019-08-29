#include "utils.hpp"

void init_logging(const std::string& log_level){
    auto g_log_level{logging::trivial::info};
    if (log_level == "debug")
        g_log_level = logging::trivial::debug;
    boost::log::core::get()->set_filter
            (
                    logging::trivial::severity >= g_log_level
            );
}

int encode_char(const char &c) {

    switch (c) {
        case 'A':
        case 'a':
            return 1;

        case 'C':
        case 'c':
            return 2;

        case 'G':
        case 'g':
            return 3;

        case 'T':
        case 't':
            return 4;

        default:
            throw std::invalid_argument("The provided character is not one of the four standard nucleotides.");
    }
}
