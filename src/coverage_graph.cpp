#include "coverage_graph.hpp"

coverage_Node::coverage_Node(std::string const seq, int const pos, int const site_ID, int const allele_ID) :
sequence(seq), pos(pos), site_ID(site_ID), allele_ID(allele_ID) {
   if (site_ID == 0 && allele_ID == 0) is_in_site = false;
   else is_in_site = true;

    coverage = std::vector<uint64_t >{sequence.length(), 0};
}

coverage_Node::coverage_Node(auto_Node const& node_in, int const site_ID, int const allele_ID) :
sequence(node_in.characters), pos(node_in.pos), site_ID(site_ID), allele_ID(allele_ID){
    if (site_ID == 0 && allele_ID == 0) is_in_site = false;
    else is_in_site = true;

    coverage = std::vector<uint64_t>{sequence.length(), 0};
}


coverage_Graph::coverage_Graph(sequence_Graph const& graph_in){
    // Throughout, we must distinguish between nodes in the constructor parameter's graph,
    // and nodes in the new graph.

    // maps the coverage graph bubble equivalent.
    std::unordered_map<std::shared_ptr<auto_Node>, std::shared_ptr<coverage_Node> > bubble_translation_map;
    auto& t_b_m = graph_in.bubble_map; // Template bubble map
    std::shared_ptr<coverage_Node> backWire;
    int cur_pos; std::string seqBuffer{""};
    int site_ID = 0;
    int allele_ID;
    std::set<std::shared_ptr<coverage_Node>> known_fixed_points;

    // Bubbles are processed children first, so when we encounter them they should be mapped already.
    for (auto& s : graph_in.bubble_map){
        allele_ID = 0;
        site_ID++;
        auto bubble_entry = std::make_shared<coverage_Node>(*s.first, site_ID, allele_ID);
        auto bubble_exit = std::make_shared<coverage_Node>(*s.second, site_ID, 0);
        cur_pos = s.first->pos + 1;
        bool skip_fixed_point{false};

        for (auto& path_seed : s.first->next){
            backWire = bubble_entry;
            allele_ID++;
           auto cur_Node = path_seed;
           while (cur_Node != s.second){
               if (cur_Node->prev.size() >= 1 && cur_Node->next.size() == 1){
                   seqBuffer += cur_Node->characters;
                   cur_Node = *(cur_Node->next.begin());
                   // If we have reached the bubble end, we need not to skip processing the sequence buffer.
                   if (cur_Node != s.second) continue;
               }

               if (seqBuffer.size() > 0){
                   auto new_Node = std::make_shared<coverage_Node>(seqBuffer, cur_pos, site_ID, allele_ID);
                   new_Node->prev.insert(backWire);
                   backWire->next.insert(new_Node);
                   backWire = new_Node;
                   seqBuffer = "";
               }

               if (t_b_m.find(cur_Node) != t_b_m.end() && cur_Node != s.second){
                   auto& translated_bubble = bubble_translation_map.at(cur_Node);
                   translated_bubble->prev.insert(backWire);
                   backWire->next.insert(translated_bubble);
                   backWire = bubble_map.at(translated_bubble);
                   cur_Node = t_b_m.at(cur_Node);
                   cur_pos = cur_Node->pos + 1;
                   // We may end at the same place as a previously processed bubble.
                   if (!skip_fixed_point && cur_Node == s.second &&
                        known_fixed_points.find(backWire) != known_fixed_points.end()){
                       skip_fixed_point = true;
                       bubble_exit = backWire;
                   }
                   continue;
               }

               // We ought to be at the bubble end point here.
               assert(cur_Node == s.second);
           }
           if (backWire != bubble_exit){
               backWire->next.insert(bubble_exit);
               bubble_exit->prev.insert(backWire);
           }
        }
        bubble_translation_map.insert(std::make_pair(s.first, bubble_entry));
        bubble_map.insert(std::make_pair(bubble_entry,bubble_exit));
        known_fixed_points.insert(bubble_exit);
    };
}
