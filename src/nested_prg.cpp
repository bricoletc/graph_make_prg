#include "nested_prg.hpp"
#include <queue>
#include <unordered_map>
#include <set>


// TODO: after joining, delete branchpoint? And how to delete its entry in map m?

// TODO: multifurcations.

prg_Node::prg_Node()
        :
        sequence(""),
        branch_point(NULL),
        multifurc(false) {
}

prg_Node::prg_Node(std::string sequence, prg_Node *branch_point)
        :
        sequence(sequence),
        branch_point(branch_point),
        multifurc(false) {
}

nested_prg::nested_prg(auto_Node *root) {
    // Set up data structures
    std::queue<auto_Node *> q;
    std::set<auto_Node *> to_visit;
    std::set<auto_Node *> visited;
    std::unordered_map<auto_Node *, prg_Node *> m;

    auto_Node *cur_Node;

    q.push(root);
    to_visit.insert(root);

    prg_Node *p_Node = new prg_Node();
    m.insert(std::make_pair(root, p_Node));

    std::string cur_Prg;

    while (!q.empty()) {
        cur_Node = q.front();
        q.pop();

        visited.insert(cur_Node);
        to_visit.erase(cur_Node);

        // Add the character of the mapped NFA node to the prg_Node.
        p_Node = m.at(cur_Node);

        // Note, that if we are in a 'join' case, p_Node->sequence ought to be empty string before this line.
        if (cur_Node -> letter != SOURCE_CHAR && cur_Node -> letter != SINK_CHAR) p_Node->sequence += cur_Node->letter;

        std::cout << "Cur_Node letter: " << cur_Node->letter << "\t at Pos: " << cur_Node->pos <<  "\t p_Node sequence: " << p_Node->sequence<<std::endl;

        // Case: join operation
        if (cur_Node->prev.size() > 1) {
            bool repushed = false;

            for (auto node : cur_Node->prev) {
                if (visited.find(node) ==
                    visited.end()) { // We need to have processed all immediate parents; if not, defer processing the node.
                    p_Node->sequence.pop_back(); // Remove the character that we added to the sequence.
                    visited.erase(cur_Node);
                    to_visit.insert(cur_Node);
                    q.push(cur_Node); // See you later
                    repushed = true;
                    break;
                }
            }

            if (repushed) continue;

            std::string join = "";
            prg_Node *branch_point;
            bool first = true; // Special treatment for first parent.
            bool multifurc = false; // Special treatment for multifurcating cases.

            prg_Node *ancestor_p_Node;

            std::cout << "JOINING: \n";
            for (auto node : cur_Node->prev) {
                ancestor_p_Node = m.at(node);

                if (p_Node->branch_point == ancestor_p_Node){
                    // Empty record
                    join += ",";
                    continue;
                }

                if (first) {
                    // Assumption: the branch_point of all the nodes we are joining, is the same- so take the first as representative.
                    branch_point = ancestor_p_Node->branch_point;
                    first = false;
                    if (branch_point->multifurc) multifurc = true;
                }

                join += ancestor_p_Node->sequence + ",";

                // Free the 'parental' prg_Node
                delete ancestor_p_Node;
                m.erase(node);

                std::cout << join << std::endl;
            }
            join.pop_back(); // Take out the last comma
            join = "[" + join + "]";
            if (branch_point != NULL && !multifurc)
                join = branch_point->sequence + join; //Only add ancestral sequence if no multifurcation occurred.
            p_Node->sequence = cur_Node->letter == SINK_CHAR ? join : join +
                                                                      p_Node->sequence; //Only add p_Node sequence if we are not at final node.
            if (branch_point != NULL) p_Node->branch_point = branch_point->branch_point; // Assign ancestor's ancestor, as new ancestor.

            std::cout << "Join output: " << join << std::endl;
            std::cout << "Prg_node seq: " << p_Node->sequence << std::endl;

        }

        // No else: join operation (backwards operation) independent of linear advance/ bifurcate/ end operation (forwards operation)
        int num_descendants = cur_Node->next.size();

        // Case: linear advance
        if (num_descendants == 1) {
            auto next_Node = *(cur_Node->next.begin());
            if (to_visit.find(next_Node) != to_visit.end())
                continue; // This node points to a multi-parent node. Nothing to do.

            // Register a future visit
            q.push(next_Node);
            to_visit.insert(next_Node);


            // Update references
            if (next_Node->prev.size() == 1) {
                // Mono-parent node: update reference (= next_Node points to what cur_Node used to)
                m.insert(std::make_pair(next_Node, p_Node));
                m.erase(cur_Node);
            } else {
                // Multi-parent node: create a blank prg_Node for it.
                // Also, no touching the reference to cur_Node in the map.
                prg_Node *new_Node = new prg_Node();
                m.insert(std::make_pair(next_Node, new_Node));
            }
        }


        //Case: bifurcate
        else if (num_descendants > 1) {
            // Do we bifurcate into a multiparent, which has already been discovered? (due to another, prior bifurcation)
            bool bifurcate_into_cognate_multiparent = false;
            prg_Node *branch_point = p_Node;

            for (auto nn : cur_Node->next) {
                if (to_visit.find(nn) != to_visit.end()) {
                    std::cout << "BIFURC_MULTIPAR" << std::endl;
                    bifurcate_into_cognate_multiparent = true; // Oh-oh...
                    p_Node->multifurc = true;
                }
            }


            prg_Node *new_Node;
            for (auto nn : cur_Node->next) {

                if (to_visit.find(nn) != to_visit.end()) continue;

                // If at least one bifurcation into a multiparent has been found, only create a prg_Node for clean nodes.
                if (bifurcate_into_cognate_multiparent) {
                    new_Node = new prg_Node(p_Node->sequence, branch_point);
                }

                else {
                    // TODO possible to avoid empty records?
                    new_Node = new prg_Node("", branch_point);
                }
                m.insert(std::make_pair(nn, new_Node));
                q.push(nn);
                to_visit.insert(nn);
            }

        }

       if (cur_Node->letter == SINK_CHAR) prg = p_Node->sequence;
    }

}
