#include "nested_prg.hpp"
#include <unordered_map>
#include <set>
#include <deque>
#include <queue>
#include <string>


prg_Node::prg_Node()
        :
        sequence(""),
        next(NULL) {
}

prg_Node::prg_Node(std::string sequence, auto_Node *next)
        :
        sequence(sequence),
        next(next) {
}

nested_prg::nested_prg(auto_Node *root) {
    map_all_bubbles(root);

    // Parse the bubbles, in dependency order
    while(!topological_order.empty()){
        auto start_point = topological_order.top();
        topological_order.pop();
        auto end_point = bubble_map.at(start_point);

        parse_bubbles(start_point,end_point);
    }

    // Linear advance to build final prg sequence

    auto_Node* cur_Node = root;
    while (cur_Node->letter != SINK_CHAR){
       if (prg_map.find(cur_Node) != prg_map.end()){
           auto p_Node = prg_map.at(cur_Node);
           prg += p_Node->sequence;
           cur_Node = p_Node->next;
       }

       else{
           if (cur_Node->letter != SOURCE_CHAR && fixed_point_map.find(cur_Node) == fixed_point_map.end()){
             prg += cur_Node->letter;
           }
           cur_Node = *(cur_Node->next.begin());
       }
    }
}


void nested_prg::map_all_bubbles(auto_Node* root){
    auto_Node *cur_Node = root;
    while (cur_Node->letter != SINK_CHAR){
        while (cur_Node->next.size() == 1) cur_Node = *(cur_Node->next.begin());

        if (cur_Node->letter == SINK_CHAR) break;
        cur_Node = map_bubbles(cur_Node);
    }
}


auto_Node *nested_prg::map_bubbles(auto_Node *start_point) {
    // Commit this start point to the topological order stack.
    topological_order.push(start_point);

    //std::deque<auto_Node*> q;
    auto cmp = [](auto_Node *p1, auto_Node *p2) { return p1->pos > p2->pos; };
    std::priority_queue<auto_Node *, std::vector<auto_Node *>, decltype(cmp)> q(cmp);
    std::set<auto_Node *> to_visit;


    for (auto nn : start_point->next) {
        q.push(nn);
        to_visit.insert(nn);
    }

    auto_Node *cur_Node;

    while (q.size() > 1) {

        cur_Node = q.top();
        q.pop();
        to_visit.erase(cur_Node);


        auto num_descendents = cur_Node->next.size();

        // Recursive call here.
        if (num_descendents > 1) {
            auto nn = map_bubbles(cur_Node);
            if (to_visit.find(nn) != to_visit.end()) continue;
            q.push(nn);
            to_visit.insert(nn);
        } else {
            auto nn = *(cur_Node->next.begin());
            if (to_visit.find(nn) != to_visit.end()) continue;
            q.push(nn);
            to_visit.insert(nn);
        }
    }

    auto fixed_point = q.top();
    // Mark the direct ancestors of the fixed point
    if (fixed_point_map.find(fixed_point) == fixed_point_map.end()){
        std::set<auto_Node*> prevs;
        for (auto pn : fixed_point->prev){
           prevs.insert(pn);
        }
        fixed_point_map.insert(std::make_pair(fixed_point, prevs));
    }

    std::cout << "Fixed point : " << fixed_point->letter << std::endl;
    bubble_map.insert(std::make_pair(start_point, fixed_point));

    return fixed_point;
}


void nested_prg::parse_bubbles(auto_Node *start_point, auto_Node *end_point) {
    std::vector<std::string> alts;
    auto ancestors = fixed_point_map.at(end_point);


    for (auto nn : start_point->next) {
        // Deal with the case where start_point has straight connection to end_point
        if (nn == end_point && ancestors.find(start_point) != ancestors.end()) ancestors.erase(start_point);

        std::string alt = "";

        while (nn != end_point) {
            if (prg_map.find(nn) != prg_map.end()) {
                auto p_Node = prg_map.at(nn);
                alt += (p_Node->sequence);
                nn = p_Node->next;
            } else {
                alt += (nn->letter);
                // Important: peel off fixed point ancestor if first time seen
                auto next = *(nn->next.begin());
                if (next == end_point && ancestors.find(nn) != ancestors.end()) ancestors.erase(nn);

                nn = next;
            }
        }

        alts.push_back(alt);
    }

    // Update ancestors
    fixed_point_map.erase(end_point);
    fixed_point_map.insert(make_pair(end_point,ancestors));

    // Make a prg sequence.
    std::string prg_Seq;
    for (auto s : alts){
       prg_Seq+= s + ",";
    }
    prg_Seq.pop_back();
    prg_Seq = "[" + prg_Seq + "]";

    // Prepend the common string pre bifurcation.
    if (start_point->letter != SOURCE_CHAR) prg_Seq = start_point->letter + prg_Seq;

    // Postpend the common letter post joining
    if (ancestors.size() == 0 && bubble_map.find(end_point) == bubble_map.end() && end_point->letter != SINK_CHAR){
      prg_Seq += end_point->letter;
    }

    std::cout << prg_Seq << std::endl;;

    // Make a prg node, and make it available.
    prg_Node* new_Node = new prg_Node(prg_Seq, end_point);
    prg_map.insert(std::make_pair(start_point,new_Node));
}
