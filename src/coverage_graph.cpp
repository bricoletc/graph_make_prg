#include "coverage_graph.hpp"

coverage_Node::coverage_Node(std::string const seq, int const pos, int const site_ID, int const allele_ID) :
sequence(seq), pos(pos), site_ID(site_ID), allele_ID(allele_ID) {
   if (site_ID == 0 && allele_ID == 0) is_in_site = false;
   else is_in_site = true;

    coverage = std::vector<uint64_t >{sequence.length(), 0};
}

coverage_Node::coverage_Node(auto_Node const& node_in, int const site_ID, int const allele_ID) :
        sequence(node_in.sequence), pos(node_in.pos), site_ID(site_ID), allele_ID(allele_ID){
    if (site_ID == 0 && allele_ID == 0) is_in_site = false;
    else is_in_site = true;

    coverage = std::vector<uint64_t>{sequence.length(), 0};
}

bool operator>(const covG_ptr &lhs, const covG_ptr &rhs) {
    if (lhs->pos == rhs->pos) {
        return lhs.get() > rhs.get();
    } else return lhs->pos > rhs->pos;
}

coverage_Graph::coverage_Graph(sequence_Graph const& graph_in){
    // Throughout, we must distinguish between nodes in the constructor parameter's graph,
    // and nodes in the new graph.

    // maps the coverage graph bubble start equivalent.
    std::unordered_map<seqG_ptr, covG_ptr > entry_translation_map;
    // maps the coverage graph bubble end equivalent.
    std::unordered_map<seqG_ptr, covG_ptr > exit_translation_map;

    covG_ptr backWire; // Points to latest previous node
    int cur_pos; std::string seqBuffer{""}; // For giving new nodes sequence
    seqG_ptr cur_Node; // For traversing the input graph.
    int site_ID = 0;
    int allele_ID;

    auto& t_b_m = graph_in.bubble_map; // constructor parameter's bubble map

    // We will copy each bubble in dependency order: child bubbles first (that is how they are ordered in the map)
    for (auto& s : graph_in.bubble_map){
        allele_ID = 0;
        site_ID++;

        // Entry and exit bubble copies.
        auto bubble_entry = std::make_shared<coverage_Node>(*s.first, site_ID, allele_ID);
        covG_ptr bubble_exit; // This one is not initialised, as it can have been made already.

        cur_pos = s.first->pos + 1;
        bool skip_fixed_point{false};
        bool deletion_bubble{false};
        std::string deletion_prefix{""};

        // Determine whether we have a direct deletion; if so we need to add some sequence
        // Meaning, we will take the last character of the bubble entry's sequence, take it out from the bubble entry,
        // And preprend it to each allele in the bubble.
        // TODO: do not do this if the bubble starts at SOURCE_CHAR
        for (auto& nn : s.first->next){
            if (nn == s.second) {
                if (!deletion_bubble){
                    deletion_bubble = true;
                    deletion_prefix = bubble_entry->sequence.back();
                    bubble_entry->sequence = bubble_entry->sequence.substr(0,bubble_entry->sequence.size() - 1);
                    cur_pos--;
                }
            }
        }

        // Traverse each allele of the bubble
        for (auto& path_seed : s.first->next){
            seqBuffer = "";
            backWire = bubble_entry;
            allele_ID++;
           cur_Node = path_seed;

           // Check for direct deletion
           // And make sure there is sequence in paths that directly join bubble start to bubble end.
           if (deletion_bubble){
              seqBuffer += deletion_prefix;
              if (cur_Node == s.second){
                  auto new_Node = std::make_shared<coverage_Node>(seqBuffer, cur_pos, site_ID, allele_ID);
                  new_Node->prev.insert(backWire);
                  backWire->next.insert(new_Node);
                  backWire = new_Node;
              }
           }

           while (cur_Node != s.second){
               if (cur_Node->prev.size() >= 1 && cur_Node->next.size() == 1){
                   // Only commit the sequence if we are not in a bubble end.
                   if (exit_translation_map.find(cur_Node) == exit_translation_map.end()){
                       seqBuffer += cur_Node->sequence;
                   }
                   cur_Node = *(cur_Node->next.begin());
                   // If we have reached the bubble end, we need to not skip processing the sequence buffer.
                   if (cur_Node != s.second) continue;
               }

               // By implication here, we are at the end point, or at a bubble start.
               // Process the sequence buffer: empty its contents into a new node.
               if (seqBuffer.size() > 0){
                   auto new_Node = std::make_shared<coverage_Node>(seqBuffer, cur_pos, site_ID, allele_ID);
                   new_Node->prev.insert(backWire);
                   backWire->next.insert(new_Node);
                   backWire = new_Node;
                   seqBuffer = "";
               }

               if (t_b_m.find(cur_Node) != t_b_m.end() && cur_Node != s.second){
                   auto translated_bubble = entry_translation_map.at(cur_Node);
                   // Avoid self link if the node is both a bubble end and a bubble start.
                   if (backWire != translated_bubble){
                       translated_bubble->prev.insert(backWire);
                       backWire->next.insert(translated_bubble);
                   }
                   backWire = bubble_map.at(translated_bubble);
                   cur_Node = t_b_m.at(cur_Node);
                   cur_pos = cur_Node->pos + 1;
                   continue;
               }

               // We must be at the bubble end point here.
               assert(cur_Node == s.second);
           }


           // We might end in a bubble that we have ended in before.
           // If we do, reuse it!
            if (exit_translation_map.find(s.second) != exit_translation_map.end()){
                bubble_exit = exit_translation_map.at(s.second);
            }
            else { // Never seen it before: make it
                // If the bubble ends where a bubble start, use that node!
                if (t_b_m.find(cur_Node) != t_b_m.end()){
                    auto translated_bubble = entry_translation_map.at(cur_Node);
                    bubble_exit = translated_bubble;
                }
                else bubble_exit = std::make_shared<coverage_Node>(*s.second, 0, 0);
                exit_translation_map.insert({s.second, bubble_exit});
                // Map the number of incidents to the bubble if first time we look at bubble end.
                int num_incidents = graph_in.fixed_point_numbers.at(s.second);
                fixed_point_numbers.insert(std::make_pair(bubble_exit, num_incidents));
            }

            // The backWire can be bubble_exit if the allele path ends by going through a bubble ending at bubble exit.
            // If that's the case, we want to avoid self linking.
            if (backWire != bubble_exit){
                backWire->next.insert(bubble_exit);
               bubble_exit->prev.insert(backWire);
           }
        }
        entry_translation_map.insert(std::make_pair(s.first, bubble_entry));
        bubble_map.insert(std::make_pair(bubble_entry,bubble_exit));
    };

    // Finally, a linear traversal to generate and wire the non-variant bits.
    seqBuffer = "";
    // Initialisation: make the root node sequence if it does not exist, or skip past it if it is a bubble start.
    cur_Node = graph_in.root;
    assert(cur_Node->sequence == SOURCE_CHAR);
    if (cur_Node->next.size() > 1){
        auto& translated_bubble = entry_translation_map.at(cur_Node); // Will throw error if not there; it ought to be.
        backWire = bubble_map.at(translated_bubble);
        cur_Node = t_b_m.at(cur_Node);
        cur_pos = cur_Node->pos + 1;
        root = translated_bubble;
    }
    else {
        cur_pos = cur_Node->pos;
        auto new_Node = std::make_shared<coverage_Node>(SOURCE_CHAR, cur_pos, 0, 0);
        backWire = new_Node;
        cur_Node = *(cur_Node->next.begin());
        root = new_Node;
    }

    while(cur_Node->sequence != SINK_CHAR){
       if (cur_Node->next.size() == 1) {
           // Only commit the sequence if we are not in a bubble end.
           if (exit_translation_map.find(cur_Node) == exit_translation_map.end()){
               seqBuffer += cur_Node->sequence;
           }
           cur_Node = *(cur_Node->next.begin());
           if (cur_Node->sequence != SINK_CHAR) continue;
       }

       if (seqBuffer.size() > 0){
          auto new_Node = std::make_shared<coverage_Node>(seqBuffer, cur_pos, 0, 0);
          seqBuffer = "";
          new_Node->prev.insert(backWire);
          backWire->next.insert(new_Node);
           backWire = new_Node;
          }

       // Case: we have a bubble
       if (cur_Node->sequence != SINK_CHAR){
           auto& translated_bubble = entry_translation_map.at(cur_Node); // Will throw error if not there; it ought to be.
           // The backWire will be the same as the bubble entry if there is a
           // node which is both a bubble end and a bubble start. In that case, no link needed.
           if (backWire != translated_bubble){
               translated_bubble->prev.insert(backWire);
               backWire->next.insert(translated_bubble);
           }
           backWire = bubble_map.at(translated_bubble);
           cur_Node = t_b_m.at(cur_Node);
           cur_pos = cur_Node->pos + 1;
       }
       }
        assert(cur_Node->sequence == SINK_CHAR);
        if (backWire->sequence != SINK_CHAR){ // The sink is not a bubble end, it does not exist yet.
            auto new_Node = std::make_shared<coverage_Node>(SINK_CHAR, cur_Node->pos, 0, 0);
            new_Node->prev.insert(backWire);
            backWire->next.insert(new_Node);
        }
    }
