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

prg_Node::prg_Node(std::string sequence, std::shared_ptr<auto_Node> next)
        :
        sequence(sequence),
        next(next) {
}

bool operator < (const incidence_fixed_point& lhs, const incidence_fixed_point& rhs){
   return lhs.pos_earliest_incident < rhs.pos_earliest_incident;
}

nested_prg::nested_prg(std::shared_ptr<auto_Node> root, int max_num_incidents)
    :
    max_num_incidents(max_num_incidents) {
    map_all_bubbles(root);
//std::cout << "Success!";
//while (!large_incidence_fixed_points.empty()){
//    auto& f = large_incidence_fixed_points.top();
//    large_incidence_fixed_points.pop();
//    std::cout << std::endl;
//    std::cout << "Pos: " << f.fixed_point->pos << std::endl;
//    std::cout << "Earliest start bubble: " << f.pos_earliest_incident << std::endl;
//    std::cout << "Num incidents: " << f.num_incidents<< std::endl;
//}
    haplotype_expand_bubbles();
    exit(0);

    // Parse the bubbles, in dependency order
    while(!topological_order.empty()){
        auto& start_point = *(topological_order.begin());
        topological_order.erase(topological_order.begin());
        auto& end_point = bubble_map.at(start_point);

        parse_bubbles(start_point,end_point);
    }

    // Linear advance to build final prg sequence
    std::shared_ptr<auto_Node> cur_Node = root;
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


void nested_prg::map_all_bubbles(std::shared_ptr<auto_Node> root){
    auto cur_Node = root;
    while (cur_Node->letter != SINK_CHAR){
        while (cur_Node->next.size() == 1) cur_Node = *(cur_Node->next.begin());

        if (cur_Node->letter == SINK_CHAR) break;
        cur_Node = map_bubbles(cur_Node);
    }

    // Populate priority queue with high incidence fixed points
    for (auto& s: fixed_point_numbers){
        if (s.second.num_incidents > max_num_incidents) large_incidence_fixed_points.insert(s.second);
    }
    fixed_point_numbers.clear(); // No longer need this.
}


std::shared_ptr<auto_Node> nested_prg::map_bubbles(std::shared_ptr<auto_Node> start_point) {
    // Return directly if we have already done this work.
    if (bubble_map.find(start_point) != bubble_map.end()) return bubble_map.at(start_point);

    // Commit this start point to the topological order stack.
    topological_order.insert(start_point);

    //std::deque<std::shared_ptr<auto_Node>> q;
    auto cmp = [](std::shared_ptr<auto_Node> p1, std::shared_ptr<auto_Node> p2) { return p1->pos > p2->pos; };
    std::priority_queue<std::shared_ptr<auto_Node> , std::vector<std::shared_ptr<auto_Node> >, decltype(cmp)> q(cmp);
    std::set<std::shared_ptr<auto_Node> > to_visit;


    for (auto nn : start_point->next) {
        q.push(nn);
        to_visit.insert(nn);
    }

    std::shared_ptr<auto_Node> cur_Node;

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
        std::set<std::shared_ptr<auto_Node>> prevs;
        for (auto pn : fixed_point->prev){
           prevs.insert(pn);
        }
        fixed_point_map.insert(std::make_pair(fixed_point, prevs));
    }

    BOOST_LOG_TRIVIAL(debug) << "Fixed point : " << fixed_point->letter << " at pos: " << fixed_point->pos;
    bubble_map.insert(std::make_pair(start_point, fixed_point));

    // Make or update the incident bubbles
    if (fixed_point_numbers.find(fixed_point) == fixed_point_numbers.end()){
        incidence_fixed_point i_fixed_point = {};
        i_fixed_point.fixed_point = fixed_point;
        i_fixed_point.earliest_incident = start_point;
        i_fixed_point.pos_earliest_incident = start_point->pos;
        i_fixed_point.haplotype_resolution = 1;
        i_fixed_point.num_incidents = 1;
       fixed_point_numbers.insert(std::make_pair(fixed_point, i_fixed_point));
    }
    else{
        auto& entry = fixed_point_numbers.at(fixed_point);
        // Update the bubble entry pointed to, if its position is before.
        if (start_point->pos < entry.pos_earliest_incident){
            entry.earliest_incident = start_point;
            entry.pos_earliest_incident = start_point->pos;
        }
        // Increment the number of incident bubbles.
        entry.num_incidents++;
    }

    return fixed_point;
}

void nested_prg::haplotype_expand_bubbles(){
    while(!large_incidence_fixed_points.empty()){
        auto& f = *(large_incidence_fixed_points.begin());
        large_incidence_fixed_points.erase(large_incidence_fixed_points.begin());
        delete_in_between(f.earliest_incident, f.fixed_point);
    }
    exit(0);
};

void nested_prg::delete_in_between(std::shared_ptr<auto_Node> start_point, std::shared_ptr<auto_Node> end_point){
   std::stack<std::shared_ptr<auto_Node>> to_visit;
   std::set<std::shared_ptr<auto_Node>> visited;
   for (auto& visitable : start_point->next){
       if (visitable != end_point) to_visit.push(visitable);
   }
   while (!to_visit.empty()){
       auto cur_node = to_visit.top();
       to_visit.pop();
       for (auto& visitable : cur_node->next){
           if (visitable != end_point && visited.find(visitable) == visited.end()){
               to_visit.push(visitable);
           }
       }
       bubble_map.erase(cur_node);
       fixed_point_map.erase(cur_node);
       topological_order.erase(cur_node);

       cur_node->next.clear();
       cur_node->prev.clear();
       visited.emplace(cur_node);
   }
   int i{0};
   for (auto& s: visited){
      i++ ;
   }
    std::cout << "Num elements visited: " << i << std::endl;
   start_point->next.clear();
   end_point->prev.clear();
}

void nested_prg::parse_bubbles(std::shared_ptr<auto_Node> start_point, std::shared_ptr<auto_Node> end_point) {
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
