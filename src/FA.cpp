#include "FA.hpp"
#include <limits>

auto_Node::auto_Node()
        :
        letter('\0'),
        fixed_Point(false),
        pos(-1){
}

auto_Node::auto_Node(char l, int pos)
        :
        letter(l),
        fixed_Point(false),
        pos(pos){
}


FA::FA(MSA &msa) {
    static std::set<char> gapping_chars = {'-', '.'};

    int pos = -1;

    root = new auto_Node(SOURCE_CHAR, pos);
    int N = msa.num_records;

    auto_Node *cur_Nodes[N];

    // Initialise to point to root.
    for (int i = 0; i < N; i++) {
        cur_Nodes[i] = root;
    }

    std::vector<char> letters;

    std::unordered_map<char, auto_Node*> new_Nodes;

    std::set<char> letters_in_column;

    while (msa.more_columns) {
        pos++;
        letters = msa.next_column();
        for (int i = 0; i < N; i++) {
            char l = letters[i];
            letters_in_column.insert(l);

            // If character indicating indel padding is found, do nothing for this letter.
            if (gapping_chars.find(l) != gapping_chars.end()) continue;
            //if (gapping_chars.find(l) != gapping_chars.end()) l = '-';

            // Create a new node; and make it accessible.
            if (new_Nodes.find(l) == new_Nodes.end()){
                auto new_Node = new auto_Node(l, pos);
                new_Nodes.insert(std::make_pair(l, new_Node));
            }

            auto_Node* retrieved_node = new_Nodes.at(l);
            auto_Node* preceding_node = cur_Nodes[i];

            // If there is no edge between these two nodes, make one
            if (preceding_node->next.find(retrieved_node) == preceding_node->next.end()){
                preceding_node->next.insert(retrieved_node);
                retrieved_node->prev.insert(preceding_node);
            };


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
    auto_Node* sink_node = new auto_Node(SINK_CHAR, sink_pos);

    for (int i = 0; i < N; i++) {
        auto_Node* cn = cur_Nodes[i];

        if (cn->next.find(sink_node) == cn->next.end()){
          cn->next.insert(sink_node);
          sink_node->prev.insert(cn);
        }

    }
    sink_node->mark_as_fixed_point();

}
