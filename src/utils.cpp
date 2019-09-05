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

std::string to_readable_string(std::vector<uint64_t> const& prg){
    std::string readable_string(prg.size(), '0');
    std::unordered_map<int,int> last_allele_indices; // Will record where to close the sites.

    int pos{-1};
    for (auto& s : prg){
        pos++;
        if (s > 4){
            if (s%2 == 1) readable_string[pos] = '[';
            else{
                readable_string[pos] = ',';
                if (last_allele_indices.find(s) != last_allele_indices.end()){
                    last_allele_indices.erase(s);
                }
                last_allele_indices.insert({s, pos});
            }
        }
        switch(s){
            case 1:
                readable_string[pos] = 'A';
                break;
            case 2:
                readable_string[pos] = 'C';
                break;
            case 3:
                readable_string[pos] = 'G';
                break;
            case 4:
                readable_string[pos] = 'T';
                break;
        }
    }

    // Close the sites
    for (auto s : last_allele_indices){
        auto pos = s.second;
        readable_string[pos] = ']';
    }
    return readable_string;
}
