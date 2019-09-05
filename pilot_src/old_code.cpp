// Convert a nested PRG string to int representation, with linear site numbering.
std::vector<uint64_t> full_prg_string_to_ints(std::string const& string_prg) {
    std::stack<int> marker_stack;
    int max_var_marker{3};
    int char_count{0};

    std::vector<uint64_t> encoded_prg(string_prg.size());
    BOOST_LOG_TRIVIAL(debug) << "Prg pre serialisation: " << string_prg;
    for (int i = 0; i<string_prg.size(); ++i){
        const auto &c = string_prg[i];

        switch(c) {
            case '[' : {
                max_var_marker += 2;
                marker_stack.push(max_var_marker);
                encoded_prg[char_count++] = max_var_marker;
                break;
            }

            case ']' : {
                assert(!marker_stack.empty());
                encoded_prg[char_count++] = marker_stack.top() + 1;
                marker_stack.pop();
                break;
            }

            case ',' : {
                assert(!marker_stack.empty());
                encoded_prg[char_count++] = marker_stack.top() + 1;
                break;
            }

            default : {
                try {
                    encoded_prg[char_count++] = encode_char(c);
                    break;
                }
                catch(std::exception& e){
                    BOOST_LOG_TRIVIAL(error) << e.what();
                    exit(1);
                }
            }
        }
    }

    BOOST_LOG_TRIVIAL(info) << "Number of sites produced: " << (max_var_marker -3 ) / 2;
    return encoded_prg;
}
