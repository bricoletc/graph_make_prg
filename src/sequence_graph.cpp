#include "sequence_graph.hpp"
#include <set>
#include <string>


bool operator < (const incidence_fixed_point& lhs, const incidence_fixed_point& rhs){
    return lhs.pos_earliest_incident < rhs.pos_earliest_incident;
}

sequence_Graph::sequence_Graph(std::shared_ptr<auto_Node> root, std::string MSA_file,
                               bool is_file, int max_num_incidents)
    :
    max_num_incidents(max_num_incidents),
    MSA_file(MSA_file),
    is_file(is_file),
    root(root){
    map_all_bubbles(root);
    if (!large_incidence_fixed_points.empty()){
        haplotype_expand_bubbles();
    }
}


void sequence_Graph::map_all_bubbles(std::shared_ptr<auto_Node> root){
    auto cur_Node = root;
    int haplotype_res = 1;
    while (cur_Node->sequence != SINK_CHAR){
        while (cur_Node->next.size() == 1) cur_Node = *(cur_Node->next.begin());

        if (cur_Node->sequence == SINK_CHAR) break;
        cur_Node = map_bubbles(cur_Node, haplotype_res);
    }
    populate_large_incidences();
}


std::shared_ptr<auto_Node> sequence_Graph::map_bubbles(std::shared_ptr<auto_Node> start_point, int haplotype_res) {

    bool site_convergence{false};
    std::set<std::shared_ptr<auto_Node> > seen_bubbles;

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


        // Look for site convergence
        if (seen_bubbles.find(cur_Node) != seen_bubbles.end()) site_convergence = true;
        for (auto& s : cur_Node->next){
            if (s->next.size() < 2) continue;
            if (seen_bubbles.find(s) != seen_bubbles.end()) {
                site_convergence = true;
            }
            else seen_bubbles.insert(s);
        }

        std::shared_ptr<auto_Node> nn;
        auto num_descendents = cur_Node->next.size();
        if (num_descendents > 1) {
            // Avoid re-doing the work if we have mapped the bubble already
            if (bubble_map.find(cur_Node) != bubble_map.end()) nn = bubble_map.at(cur_Node);
            // Recursive call here.
            else nn = map_bubbles(cur_Node, haplotype_res);
        }
        else if (num_descendents == 0) continue;
        else nn = *(cur_Node->next.begin());

        // Avoid processing more than once
        if (to_visit.find(nn) != to_visit.end()) continue;
        q.push(nn);
        to_visit.insert(nn);
    }

    auto fixed_point = q.top();
    // Remove site convergence under a specific edge case
    if (seen_bubbles.size() == 1 && *(seen_bubbles.begin()) == fixed_point) site_convergence = false;

    if (fixed_point_numbers.find(fixed_point) == fixed_point_numbers.end()){
        fixed_point_numbers.insert(std::make_pair(fixed_point, 1));
    } else fixed_point_numbers.at(fixed_point)++;

    BOOST_LOG_TRIVIAL(debug)<< "Fixed point : " << fixed_point->sequence << " at pos: " << fixed_point->pos;
    bubble_map.insert(std::make_pair(start_point, fixed_point));

    // Make or update the incident bubbles
    if (fixed_point_incidence_map.find(fixed_point) == fixed_point_incidence_map.end()){
        incidence_fixed_point i_fixed_point = {};
        i_fixed_point.fixed_point = fixed_point;
        i_fixed_point.earliest_incident = start_point;
        i_fixed_point.pos_earliest_incident = start_point->pos;
        i_fixed_point.haplotype_resolution = haplotype_res;
        site_convergence ? i_fixed_point.num_incidents = max_num_incidents + 1 : i_fixed_point.num_incidents = 1;
       fixed_point_incidence_map.insert(std::make_pair(fixed_point, i_fixed_point));
    }
    else{
        auto& entry = fixed_point_incidence_map.at(fixed_point);
        // Update the bubble entry pointed to, if its position is before.
        if (start_point->pos < entry.pos_earliest_incident){
            entry.earliest_incident = start_point;
            entry.pos_earliest_incident = start_point->pos;
        }
        // Increment the number of incident bubbles.
        site_convergence ? entry.num_incidents = max_num_incidents + 1 : entry.num_incidents++;
    }

    return fixed_point;
}

void sequence_Graph::haplotype_expand_bubbles(){
    MSA msa(MSA_file, is_file);
    while(!large_incidence_fixed_points.empty()){
        auto large_incidence = *(large_incidence_fixed_points.begin());
        fixed_point_incidence_map.erase(large_incidence.fixed_point); // To avoid infinite cycling in this loop

        // Is there still room for expansion?
        auto not_fully_expanded = large_incidence.haplotype_resolution <
            (large_incidence.fixed_point->pos - large_incidence.earliest_incident->pos - 1);
        // This is an edge case due to sink pos being max_int and source pos being -1 or 0, so above test fails due to
        // max_int overflow.
        auto edge_condition = large_incidence.fixed_point->sequence == SINK_CHAR &&
                large_incidence.pos_earliest_incident <= 0;

        if (not_fully_expanded || edge_condition){
            delete_in_between(large_incidence.earliest_incident, large_incidence.fixed_point);
            int haplotype_res = rebuild_in_between(msa, large_incidence);
            bubble_map.erase(large_incidence.earliest_incident); //Otherwise, will not map the bubble
            fixed_point_numbers.erase(large_incidence.fixed_point);
            map_bubbles(large_incidence.earliest_incident, haplotype_res);
        }
        else {
            BOOST_LOG_TRIVIAL(warning) << "Could not expand bubble: " <<
                large_incidence.earliest_incident << "beyond its sequence width";
        }

        BOOST_LOG_TRIVIAL(debug) << "Resolution: " << large_incidence.haplotype_resolution << "\t";
        large_incidence_fixed_points.erase(large_incidence_fixed_points.begin());

        if (large_incidence_fixed_points.empty()) populate_large_incidences();
    }
};

void sequence_Graph::delete_in_between(std::shared_ptr<auto_Node> start_point, std::shared_ptr<auto_Node> end_point){
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

       // Remove pointers where they might be
       bubble_map.erase(cur_node);
       fixed_point_numbers.erase(cur_node);
       if (fixed_point_incidence_map.find(cur_node) != fixed_point_incidence_map.end()){
           auto large_incidence = fixed_point_incidence_map.at(cur_node);
           fixed_point_incidence_map.erase(cur_node);
           large_incidence_fixed_points.erase(large_incidence);
       }

       cur_node->next.clear();
       cur_node->prev.clear();
       visited.emplace(cur_node);
   }
   int i{0};
   for (auto& s: visited){
      i++ ;
   }
    BOOST_LOG_TRIVIAL(debug) << "Num elements visited: " << i << std::endl;
   start_point->next.clear();
   end_point->prev.clear();
}

int sequence_Graph::rebuild_in_between(MSA& msa, incidence_fixed_point& i){
    // Double the haplotype resolution
    i.haplotype_resolution *= 2;
  FA fa(msa, i.earliest_incident, i.fixed_point, i.haplotype_resolution);
  return i.haplotype_resolution;
}

void sequence_Graph::populate_large_incidences(){
    for (auto& s: fixed_point_incidence_map){
        if (s.second.num_incidents > max_num_incidents) large_incidence_fixed_points.insert(s.second);
    }
}


