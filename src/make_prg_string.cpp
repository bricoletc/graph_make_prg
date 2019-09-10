#include "make_prg_string.hpp"

template<typename T>
prg_Node<T>::prg_Node()
        :
        sequence({0}),
        next(NULL) {
}

template<typename T>
prg_Node<T>::prg_Node(std::vector<uint64_t> sequence, std::shared_ptr<T> next)
        :
        sequence(sequence),
        next(next) {
}

template<typename graph_T, typename node_T>
stringified_PRG<graph_T, node_T>::stringified_PRG(graph_T& graph_in){
    bubble_map = graph_in.bubble_map;
    fixed_point_numbers = graph_in.fixed_point_numbers;
    int site_ID = 5;

    // Parse the bubbles, in dependency order
    for (auto& pair : bubble_map){
        auto start_point = pair.first;
        auto end_point = pair.second;
        parse_bubble(start_point, end_point, site_ID);
        site_ID +=2;
    }

    // Linear advance to build final prg sequence
    std::shared_ptr<node_T> cur_Node = graph_in.root;
    while (cur_Node->sequence != SINK_CHAR){
        if (prg_map.find(cur_Node) != prg_map.end()){
            auto p_Node = prg_map.at(cur_Node);
            prg.insert(prg.end(), p_Node->sequence.begin(), p_Node->sequence.end());
            cur_Node = p_Node->next;
        }

        else{
            if (cur_Node->sequence != SOURCE_CHAR && fixed_point_numbers.find(cur_Node) == fixed_point_numbers.end()){
                auto new_sequence = prg_string_to_ints(cur_Node->sequence);
                prg.insert(prg.end(), new_sequence.begin(), new_sequence.end());
            }
            cur_Node = *(cur_Node->next.begin());
        }
    }
};


template<typename graph_T, typename node_T>
void stringified_PRG<graph_T,node_T>::parse_bubble(std::shared_ptr<node_T> start_point,
        std::shared_ptr<node_T> end_point, int site_ID){
        std::vector<std::vector<uint64_t>> alts;
        auto& num_bubbles_to_process = fixed_point_numbers.at(end_point);
        bool direct_deletion = false;
        std::set<std::shared_ptr<node_T>> used_sites;

        for (auto nn : start_point->next) {
            std::vector<uint64_t> alt;

            while (nn != end_point) {
                if (prg_map.find(nn) != prg_map.end()) {
                    auto p_Node = prg_map.at(nn);
                    // The below assertion is CRUCIAL: ideally, we want each bubble in the graph to be
                    // associated with one, and one only, site ID.
                    // This amounts to not using the same site more than ONCE in the PRG string.
                    assert(used_sites.find(nn) == used_sites.end());
                    used_sites.insert(nn);
                    alt.insert(alt.end(),p_Node->sequence.begin(),p_Node->sequence.end());
                    nn = p_Node->next;
                } else {
                    try{ // Has this sequence been committed to PRG string already?
                        if (fixed_point_numbers.at(nn) == 0 ) {;} // Yes: Do nothing
                    }
                    catch(const std::out_of_range &e) { //No: commit the sequence
                        auto new_seq = prg_string_to_ints(nn->sequence);
                        alt.insert(alt.end(),new_seq.begin(),new_seq.end());
                    }
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

        // Prepend the common string pre bifurcation.
        if (start_point->sequence != SOURCE_CHAR && !direct_deletion){
            auto new_seq = prg_string_to_ints(start_point->sequence);
            prg_Seq.insert(prg_Seq.end(),new_seq.begin(), new_seq.end());
        }

        // Site entry marker
        prg_Seq.emplace_back(site_ID);

        // Sort the alleles, so that unit testing is simplified.
        std::sort(alts.begin(), alts.end());
        // Add the alleles
        for (auto s : alts){
            // Prepend the common string to each allele if direct deletion
            if (direct_deletion){
                auto new_seq = prg_string_to_ints(start_point->sequence);
                s.insert(s.begin(),new_seq.begin(), new_seq.end());
            }
            prg_Seq.insert(prg_Seq.end(),s.begin(), s.end());
            prg_Seq.emplace_back(site_ID + 1);
        }

        // Postpend the common sequence post joining
        // It is OK to postpend at the last incident bubble, because the bubbles are ordered topologically (innermost first)
        if (num_bubbles_to_process == 0 && bubble_map.find(end_point) == bubble_map.end() && end_point->sequence != SINK_CHAR){
            auto new_seq = prg_string_to_ints(end_point->sequence);
            prg_Seq.insert(prg_Seq.end(), new_seq.begin(), new_seq.end());
        }

        // BOOST_LOG_TRIVIAL(debug) << prg_Seq ;

        // Make a prg node, and make it available.
        prg_Node<node_T>* new_Node = new prg_Node<node_T>(prg_Seq, end_point);
        prg_map.insert(std::make_pair(start_point,new_Node));
    }

std::vector<uint64_t> prg_string_to_ints(std::string const& string_prg){
    std::vector<uint64_t> encoded_prg(string_prg.size());
    int char_count{0};
    for (auto& c : string_prg){
        try {
            encoded_prg[char_count++] = encode_char(c);
        }
        catch(std::exception& e){
            BOOST_LOG_TRIVIAL(error) << e.what();
            exit(1);
        }
    }
    return encoded_prg;
}

// Define the template specialisations that can be, for linkage to occur
template class stringified_PRG<sequence_Graph, auto_Node>;
template class stringified_PRG<coverage_Graph, coverage_Node>;

