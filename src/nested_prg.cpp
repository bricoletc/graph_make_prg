#include "nested_prg.hpp"
#include <unordered_map>
#include <set>
#include <string>
#include <stack>


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
    std::cout << "Success!";
    for (auto f: fixed_point_numbers){
        if (f.second.second <= 3 or f.first->pos > 290) continue;
        std::cout << std::endl;
        std::cout << "Pos: " << f.first->pos << std::endl;
        std::cout << "Earliest start bubble: " << f.second.first->pos << std::endl;
        std::cout << "Num incidents: " << f.second.second << std::endl;
    }
    exit(0);

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

    serialise_prg();
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
    // Return directly if we have already done this work.
    if (bubble_map.find(start_point) != bubble_map.end()) return bubble_map.at(start_point);

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

    BOOST_LOG_TRIVIAL(debug) << "Fixed point : " << fixed_point->letter << " at pos: " << fixed_point->pos;
    bubble_map.insert(std::make_pair(start_point, fixed_point));

    // Update the incident bubbles
    if (fixed_point_numbers.find(fixed_point) == fixed_point_numbers.end()){
       fixed_point_numbers.insert(std::make_pair(fixed_point, std::make_pair(start_point, 1)));
    }
    else{
        auto& entry = fixed_point_numbers.at(fixed_point);
        auto& earliest_incident_bubble = entry.first;
        // Update the bubble entry pointed to, if its position is before.
        if (start_point->pos < earliest_incident_bubble->pos) earliest_incident_bubble = start_point;
        // Increment the number of incident bubbles.
        entry.second++;
    }

    return fixed_point;
}


void nested_prg::parse_bubbles(auto_Node *start_point, auto_Node *end_point) {
    std::vector<std::string> alts;
    auto ancestors = fixed_point_map.at(end_point);


    for (auto nn : start_point->next) {
        std::string alt = "";
        // Record the previous node, to wipe it out as an ancestor at the end of the loop.
        auto previous = start_point;

        while (nn != end_point) {
            previous = nn;
            if (prg_map.find(nn) != prg_map.end()) {
                auto p_Node = prg_map.at(nn);
                alt += (p_Node->sequence);
                nn = p_Node->next;
            } else {
                try{ // Has this char been committed to PRG string already?
                    if (fixed_point_map.at(nn).size() == 0 ) {;} // Yes: Do nothing
                }
                catch(const std::out_of_range &e) {alt += (nn->letter);} //No: commit the char
                nn = *(nn->next.begin());
            }
        }
        // Important!!: peel off fixed point direct ancestor first time seen.
        if (ancestors.find(previous) != ancestors.end()) ancestors.erase(previous);

        alts.push_back(alt);
    }

    // Update ancestors
    fixed_point_map.erase(end_point);
    fixed_point_map.insert(make_pair(end_point,ancestors));

    // Make a prg sequence.
    std::string prg_Seq;

    // Make sure we sort the alleles, so that unit testing is simplified.
    std::sort(alts.begin(), alts.end());
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

    BOOST_LOG_TRIVIAL(debug) << prg_Seq ;

    // Make a prg node, and make it available.
    prg_Node* new_Node = new prg_Node(prg_Seq, end_point);
    prg_map.insert(std::make_pair(start_point,new_Node));
}

void nested_prg::serialise_prg() {
   std::stack<int> marker_stack;
   int max_var_marker{3};

   std::string serialised_prg = "";
   BOOST_LOG_TRIVIAL(debug) << "Prg pre serialisation: " << prg;
   for (int i = 0; i<prg.size(); ++i){
       const auto &c = prg[i];

      switch(c) {
          case '[' : {
              max_var_marker += 2;
              marker_stack.push(max_var_marker);
              serialised_prg += " " + std::to_string(max_var_marker) + " ";
              break;
          }

          case ']' : {
              assert(!marker_stack.empty());
              serialised_prg += std::to_string(marker_stack.top() + 1) + "$ ";
              marker_stack.pop();
              break;
          }

          case ',' : {
              assert(!marker_stack.empty());
              serialised_prg += std::to_string(marker_stack.top() + 1);
              break;
          }

          default : serialised_prg += c;
                    break;
      }
   }

   this->serialised_prg = serialised_prg;
   BOOST_LOG_TRIVIAL(info) << "Number of sites produced: " << (max_var_marker -3 ) / 2;
}
