#include "nested_prg.hpp"
#include <queue>
#include <unordered_map>
#include <set>

prg_Node::prg_Node()
:
sequence(""),
branch_point(NULL){
}

prg_Node::prg_Node(std::string sequence, prg_Node *branch_point)
:
sequence(sequence),
branch_point(branch_point){
}

nested_prg::nested_prg(auto_Node *root) {
    // Set up data structures
    std::queue<auto_Node*> q;
    std::set<auto_Node*> to_visit;
    std::set<auto_Node*> visited;
    std::unordered_map<auto_Node*,prg_Node*> m;

    auto_Node* cur_Node;

    q.push(root);
    to_visit.insert(root);

    prg_Node* p_Node = new prg_Node();
    m.insert(std::make_pair(root, p_Node));

    while(!q.empty()){
        cur_Node = q.front();
        q.pop();

        visited.insert(cur_Node);
        to_visit.erase(cur_Node);

        // Add the character of the mapped NFA node to the prg_Node.
        p_Node = m.at(cur_Node);
        p_Node->sequence += cur_Node->letter; // Note, that if we are in a 'join' case, p_Node->sequence ought to be empty string before this line.

        // Case: join operation
        if (cur_Node->prev.size() > 1){
           for (auto node : cur_Node->prev){
               if (visited.find(node) == visited.end()){ // We need to have processed all immediate parents; if not, defer processing the node.
                   p_Node->sequence.pop_back(); // Remove the character that we added to the sequence.
                   q.push(cur_Node); // See you later
                   continue;
               }
           }

           std::string join = "";
           prg_Node* branch_point;
           bool first = true; // Special treatment for first parent.

           prg_Node* ancestor_p_Node;

           for (auto node : cur_Node->prev){
               ancestor_p_Node = m.at(cur_Node);
               join += ancestor_p_Node->sequence + ",";

               if (first){
                   // Assumption: the branch_point of all the nodes we are joining, is the same- so take the first as representative.
                 branch_point = ancestor_p_Node->branch_point;
                 join += ancestor_p_Node->sequence; // First allele; no preceding comma
                 first = false;
               }
               else join += "," + ancestor_p_Node->sequence;

               // Free the 'parental' prg_Node
               delete ancestor_p_Node;
               m.erase(node);
           }
           join = "[" + join + "]";
           join = branch_point->sequence + join;
           p_Node->sequence = cur_Node->letter == SINK_CHAR ? join : join + p_Node->sequence; //Only add p_Node sequence if we are not at final node.
           p_Node->branch_point = branch_point->branch_point; // Assign ancestor's ancestor, as new ancestor.
        }

        // No else: join operation (backwards operation) independent of linear advance/ bifurcate/ end operation (forwards operation)
        int num_descendants = cur_Node->next.size();

        // Case: linear advance
        if (num_descendants == 1){
            auto next_Node = *(cur_Node->next.begin());
            if (to_visit.find(next_Node) != to_visit.end()) continue; // This node points to a multi-parent node. Nothing to do.

            // Register a future visit
            q.push(next_Node);
            to_visit.insert(next_Node);


            // Update references
            if (next_Node->prev.size() == 1){
                // Mono-parent node: update reference (= next_Node points to what cur_Node used to)
                m.insert(std::make_pair(next_Node, p_Node));
                m.erase(cur_Node);
            }
            else {
                // Multi-parent node: create a blank prg_Node for it.
                // Also, don't touch the reference to cur_Node in the map.
                prg_Node* new_Node = new prg_Node();
                m.insert(std::make_pair(next_Node, new_Node));
            }
        }

        //Case: bifurcate
        else if (num_descendants > 1){
            bool bifurcate_into_multiparent = false;
            int num_clean = 0;
            prg_Node* branch_point = NULL;

            for (auto nn : cur_Node->next){
                if (to_visit.find(nn) != to_visit.end()){
                    bifurcate_into_multiparent = true; // Oh-oh...
                }
                else num_clean++;
            }

            if (num_clean > 1) branch_point = p_Node;

            for (auto nn : cur_Node->next){
                // If at least one bifurcation into a multiparent has been found, only create a prg_Node for clean nodes.
               if (bifurcate_into_multiparent && to_visit.find(nn) == to_visit.end()){
                  prg_Node* new_Node = new prg_Node(p_Node->sequence,branch_point);
                  m.insert(std::make_pair(nn, new_Node));
               }

               // The branch point should not be null here, as there will be more >1 clean nodes (because there are 0 unclean ones).
               else if (!bifurcate_into_multiparent){
                   prg_Node* new_Node = new prg_Node("", branch_point);
                   m.insert(std::make_pair(nn, new_Node));
               }
            }
        }

    }
}
