#include "FA.hpp"
#include <limits>

auto_Node::auto_Node()
        :
        characters(""),
        fixed_Point(false),
        pos(-1){
}

auto_Node::auto_Node(std::string l, int pos)
        :
        characters(l),
        fixed_Point(false),
        pos(pos){
}

bool operator > (const std::shared_ptr<auto_Node>& lhs, const std::shared_ptr<auto_Node>& rhs){
    if (lhs->pos == rhs->pos){
    return lhs.get() > rhs.get();
    }
    else return lhs->pos > rhs->pos;
}

FA::FA(MSA &msa) {
    int pos = -1;

    root = std::make_shared<auto_Node>(SOURCE_CHAR, pos);
    int N = msa.num_records;

    std::shared_ptr<auto_Node> cur_Nodes[N];

    // Initialise to point to root.
    for (int i = 0; i < N; i++) {
        cur_Nodes[i] = root;
    }

    std::vector<char> letters;

    std::unordered_map<char, std::shared_ptr<auto_Node>> new_Nodes;

    std::set<char> letters_in_column;

    while (msa.more_columns) {
        pos++;
        letters = msa.next_column();
        for (int i = 0; i < N; i++) {
            char l = letters[i];
            letters_in_column.insert(l);

            // If character indicating indel padding is found, do nothing for this characters.
            if (gapping_chars.find(l) != gapping_chars.end()) continue;
            //if (gapping_chars.find(l) != gapping_chars.end()) l = '-';

            // Create a new node; and make it accessible.
            if (new_Nodes.find(l) == new_Nodes.end()){
                auto new_Node = std::make_shared<auto_Node>(std::string(1,l), pos);
                new_Nodes.insert(std::make_pair(l, new_Node));
            }

            auto retrieved_node = new_Nodes.at(l);
            auto preceding_node = cur_Nodes[i];

            // If there is no edge between these two nodes, make them
            preceding_node->next.insert(retrieved_node);
            retrieved_node->prev.insert(preceding_node);

            // Make retrieved node accessible to next column.
            cur_Nodes[i] = retrieved_node;
        }

        if (letters_in_column.size() == 1){
            cur_Nodes[0]->mark_as_fixed_point();
        }

        new_Nodes.clear();
        letters_in_column.clear();
    }

    // Final link to a fixed_point auto_Node
    int sink_pos = std::numeric_limits<int>::max(); // We need sink_pos to be processed latest in the priority queue in nested prg construction.
    auto sink_node = std::make_shared<auto_Node>(SINK_CHAR, sink_pos);

    for (int i = 0; i < N; i++) {
        auto cn = cur_Nodes[i];

        if (cn->next.find(sink_node) == cn->next.end()){
          cn->next.insert(sink_node);
          sink_node->prev.insert(cn);
        }

    }
    sink_node->mark_as_fixed_point();

}

FA::FA(MSA& msa, std::shared_ptr<auto_Node> start_point,
        std::shared_ptr<auto_Node> end_point, int haplotypic_resolution){
    /**
     * Initialisation
     */
    assert(start_point->characters.length() == 1);

    // Filter for the positions of the records of interest.
    msa.reposition(start_point->pos); //pos is the next position fetched by call to `next_column()`
    std::vector<int> records_of_interest;
    int N = msa.num_records;
    auto letters = msa.next_column();
    char l;
    for (int i = 0; i < N; i++){
        l = letters[i];
        if (start_point->characters == SOURCE_CHAR ||
            std::string(1,l) == start_point->characters) records_of_interest.push_back(i);
    }
    N = static_cast<int>(records_of_interest.size());


    int pos = start_point->pos + 1; // Now, we want to start by getting the first column after the start point.
    msa.reposition(pos); //pos is the next position fetched by call to `next_column()`

    std::shared_ptr<auto_Node> cur_Nodes[N];

    // Initialise to point to start point.
    for (int i = 0; i < N; i++) {
        cur_Nodes[i] = start_point;
    }

    std::string buffer[N]; // Holds the growing sequences
    for (int i = 0; i < N; i++) {
        buffer[i] = "";
    }

    std::unordered_map<std::string, std::shared_ptr<auto_Node>> new_Nodes;

    /**
     * Iteration: build nodes & edges
     */
     bool gapping;
     int node_pos;

     int end_point_pos;
     if (end_point->characters == SINK_CHAR) end_point_pos = msa.record_size - 1;
     else end_point_pos = end_point->pos;
    while (pos != end_point_pos){
       auto letters = msa.next_column();

        for (int i = 0; i < N; i++) {
            gapping = false;
           l =  letters[records_of_interest[i]];
           if (gapping_chars.find(l) != gapping_chars.end()){
               if (buffer[i].length() == 0) continue;
               gapping = true;
               node_pos = pos - 1;
           }
           else {
               buffer[i] += l;
               node_pos = pos;
           }

           // Near the end point, collect sequence of size up to 2*haplotypic_resolution - 1 characters.
           if (!gapping && end_point->pos - pos < haplotypic_resolution*2 - 1) continue;

           // Case: we want a node.
           if (gapping || (buffer[i].length() == haplotypic_resolution)){
              if (new_Nodes.find(buffer[i]) == new_Nodes.end()){
                  auto new_Node = std::make_shared<auto_Node>(buffer[i],
                          node_pos - buffer[i].length() + 1);
                  new_Nodes.insert(std::make_pair(buffer[i], new_Node));
              }
              auto& retrieved_Node = new_Nodes.at(buffer[i]);
              auto& cur_Node = cur_Nodes[i];

              cur_Node->next.insert(retrieved_Node);
              retrieved_Node->prev.insert(cur_Node);

              cur_Nodes[i] = retrieved_Node;
              buffer[i] = "";
           }
        }

        pos++;
        new_Nodes.clear();
    }

    // End condition: link to the end_point
    std::set<std::shared_ptr<auto_Node>> final_Nodes;
    for (int i = 0; i < N; i++) {
        auto cn = cur_Nodes[i];
        if (buffer[i].length() == 0){
            final_Nodes.insert(cn);
        }
        else {
            // Case: the buffer contains non-gapping sequence- if there was a gap prior, it got committed.
            // So its start position is then end_point position less the sequence length.
            if (new_Nodes.find(buffer[i]) == new_Nodes.end()){
                auto new_Node = std::make_shared<auto_Node>(buffer[i],
                        pos - buffer[i].length());
                new_Nodes.insert(std::make_pair(buffer[i], new_Node));
                final_Nodes.insert(new_Node);
            }
            auto& new_Node = new_Nodes.at(buffer[i]);
            cn->next.insert(new_Node);
            new_Node->prev.insert(cn);
        }
    }
    for (auto& final_Node : final_Nodes){
        final_Node->next.insert(end_point);
        end_point->prev.insert(final_Node);
    }

}
