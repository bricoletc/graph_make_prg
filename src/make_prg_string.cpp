#include "make_prg_string.hpp"

prg_Node::prg_Node()
        :
        sequence({0}),
        next(NULL) {
}

prg_Node::prg_Node(std::vector<uint64_t> sequence, std::shared_ptr<auto_Node> next)
        :
        sequence(sequence),
        next(next) {
}

template<typename T>
stringified_PRG<T>::stringified_PRG(T& graph_in){
    bubble_map = &graph_in.bubble_map;
    fixed_point_numbers = &graph_in.fixed_point_numbers;

    // Parse the bubbles, in dependency order
    for (auto& pair : bubble_map){
        auto& start_point = pair.first;
        auto& end_point = pair.second;
        parse_bubble(start_point, end_point);
    }

    // Linear advance to build final prg sequence
    std::shared_ptr<T> cur_Node = graph_in.root;
    while (cur_Node->sequence != SINK_CHAR){
        if (prg_map.find(cur_Node) != prg_map.end()){
            auto p_Node = prg_map.at(cur_Node);
            prg += p_Node->sequence;
            cur_Node = p_Node->next;
        }

        else{
            if (cur_Node->sequence != SOURCE_CHAR && fixed_point_numbers.find(cur_Node) == fixed_point_numbers.end()){
                prg += cur_Node->sequence;
            }
            cur_Node = *(cur_Node->next.begin());
        }
    }
};


template<typename T>
void stringified_PRG<T>::parse_bubble(std::shared_ptr<T> start_point, std::shared_ptr<T> end_point, int site_ID){
        std::vector<std::vector<uint64_t>> alts;
        auto& num_bubbles_to_process = fixed_point_numbers.at(end_point);
        bool direct_deletion = false;
        std::set<std::shared_ptr<T>> used_sites;

        for (auto nn : start_point->next) {
            std::vector<uint64_t> alt;

            while (nn != end_point) {
                if (prg_map.find(nn) != prg_map.end()) {
                    auto p_Node = prg_map.at(nn);
                    // The below assertion is a goal we need to achieve for sites to be unambiguously identified.
                    assert(used_sites.find(nn) == used_sites.end());
                    used_sites.insert(nn);
                    alt.emplace_back(p_Node->sequence);
                    nn = p_Node->next;
                } else {
                    try{ // Has this char been committed to PRG string already?
                        if (fixed_point_numbers.at(nn) == 0 ) {;} // Yes: Do nothing
                    }
                    catch(const std::out_of_range &e) {alt.emplace_back(nn->sequence);} //No: commit the char
                    nn = *(nn->next.begin());
                }
            }
            // This ought never to evaluate to true for `coverage_Graph`s
            if (alt.size() == 0)  direct_deletion = true;
            alts.push_back(alt);
        }

        num_bubbles_to_process -= 1;
        // Make a prg sequence.
        std::vector<uint64_t> prg_Seq;
        prg_Seq.emplace_back(site_ID);

        // Make sure we sort the alleles, so that unit testing is simplified.
        std::sort(alts.begin(), alts.end());
        for (auto s : alts){
            // Prepend the common string to each allele
            if (direct_deletion) s = s.insert(s.begin(),start_point->sequence.begin(), start_point->sequence.end());
            prg_Seq.emplace_back(site_ID + 1);
        }
        prg_Seq.pop_back();
        prg_Seq.emplace_back(site_ID);

        // Prepend the common string pre bifurcation.
        if (start_point->sequence != SOURCE_CHAR && !direct_deletion) prg_Seq = start_point->sequence + prg_Seq;

        // Postpend the common sequence post joining
        // It is OK to postpend at the last incident bubble, because the bubbles are ordered topologically (innermost first)
        if (num_bubbles_to_process == 0 && bubble_map.find(end_point) == bubble_map.end() && end_point->sequence != SINK_CHAR){
            prg_Seq += end_point->sequence;
        }

        BOOST_LOG_TRIVIAL(debug) << prg_Seq ;

        // Make a prg node, and make it available.
        prg_Node* new_Node = new prg_Node(prg_Seq, end_point);
        prg_map.insert(std::make_pair(start_point,new_Node));
    }
};


void binary_encode() {
    std::stack<int> marker_stack;
    int max_var_marker{3};
    int char_count{0};

    sdsl::int_vector<> encoded_prg(this->prg.length(), 0, 32);
    BOOST_LOG_TRIVIAL(debug) << "Prg pre serialisation: " << prg;
    for (int i = 0; i<prg.size(); ++i){
        const auto &c = prg[i];

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

    this->encoded_prg = encoded_prg;
    BOOST_LOG_TRIVIAL(info) << "Number of sites produced: " << (max_var_marker -3 ) / 2;
}

