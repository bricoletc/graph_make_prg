#include "nested_prg.hpp"
#include <set>
#include <string>


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

nested_prg::nested_prg(std::shared_ptr<auto_Node> root, std::string MSA_file,
        bool is_file, int max_num_incidents)
    :
    max_num_incidents(max_num_incidents),
    MSA_file(MSA_file),
    is_file(is_file){
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
    if (!large_incidence_fixed_points.empty()){
        haplotype_expand_bubbles();
    }

    // Parse the bubbles, in dependency order
    while(!topological_order.empty()){
        auto& start_point = *(topological_order.begin());
        topological_order.erase(topological_order.begin());
        auto& end_point = bubble_map.at(start_point);

        parse_bubbles(start_point,end_point);
    }

    // Linear advance to build final prg sequence
    std::shared_ptr<auto_Node> cur_Node = root;
    while (cur_Node->characters != SINK_CHAR){
       if (prg_map.find(cur_Node) != prg_map.end()){
           auto p_Node = prg_map.at(cur_Node);
           prg += p_Node->sequence;
           cur_Node = p_Node->next;
       }

       else{
           if (cur_Node->characters != SOURCE_CHAR && fixed_point_numbers.find(cur_Node) == fixed_point_numbers.end()){
             prg += cur_Node->characters;
           }
           cur_Node = *(cur_Node->next.begin());
       }
    }

    add_site_numbers_and_binary_encode();
}


void nested_prg::map_all_bubbles(std::shared_ptr<auto_Node> root){
    auto cur_Node = root;
    int haplotype_res = 1;
    while (cur_Node->characters != SINK_CHAR){
        while (cur_Node->next.size() == 1) cur_Node = *(cur_Node->next.begin());

        if (cur_Node->characters == SINK_CHAR) break;
        cur_Node = map_bubbles(cur_Node, haplotype_res);
    }
    populate_large_incidences();
}


std::shared_ptr<auto_Node> nested_prg::map_bubbles(std::shared_ptr<auto_Node> start_point, int haplotype_res) {
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
            auto nn = map_bubbles(cur_Node, haplotype_res);
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
    if (fixed_point_numbers.find(fixed_point) == fixed_point_numbers.end()){
        fixed_point_numbers.insert(std::make_pair(fixed_point, 1));
    } else fixed_point_numbers.at(fixed_point)++;

    BOOST_LOG_TRIVIAL(debug)<< "Fixed point : " << fixed_point->characters << " at pos: " << fixed_point->pos;
    bubble_map.insert(std::make_pair(start_point, fixed_point));

    // Make or update the incident bubbles
    if (fixed_point_incidence_map.find(fixed_point) == fixed_point_incidence_map.end()){
        incidence_fixed_point i_fixed_point = {};
        i_fixed_point.fixed_point = fixed_point;
        i_fixed_point.earliest_incident = start_point;
        i_fixed_point.pos_earliest_incident = start_point->pos;
        i_fixed_point.haplotype_resolution = haplotype_res;
        i_fixed_point.num_incidents = 1;
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
        entry.num_incidents++;
    }

    return fixed_point;
}

void nested_prg::haplotype_expand_bubbles(){
    MSA msa(MSA_file, is_file);
    while(!large_incidence_fixed_points.empty()){
        auto large_incidence = *(large_incidence_fixed_points.begin());
        fixed_point_incidence_map.erase(large_incidence.fixed_point); // To avoid infinite cycling in this loop

        // Is there still room for expansion?
        auto not_fully_expanded = large_incidence.haplotype_resolution <
            (large_incidence.fixed_point->pos - large_incidence.earliest_incident->pos + 1);
        // This is an edge case due to sink pos being max_int and source pos being -1, so above test fails.
        auto edge_condition = large_incidence.fixed_point->characters == SINK_CHAR &&
                large_incidence.earliest_incident->characters == SOURCE_CHAR;

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

       // Remove pointers where they might be
       bubble_map.erase(cur_node);
       fixed_point_numbers.erase(cur_node);
       topological_order.erase(cur_node);
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

int nested_prg::rebuild_in_between(MSA& msa, incidence_fixed_point& i){
    // Double the haplotype resolution
    i.haplotype_resolution *= 2;
  FA fa(msa, i.earliest_incident, i.fixed_point, i.haplotype_resolution);
  return i.haplotype_resolution;
}

void nested_prg::populate_large_incidences(){
    for (auto& s: fixed_point_incidence_map){
        if (s.second.num_incidents > max_num_incidents) large_incidence_fixed_points.insert(s.second);
    }
}

void nested_prg::parse_bubbles(std::shared_ptr<auto_Node> start_point, std::shared_ptr<auto_Node> end_point) {
    std::vector<std::string> alts;
    auto& num_bubbles_to_process = fixed_point_numbers.at(end_point);
    bool direct_deletion = false;
    std::set<std::shared_ptr<auto_Node>> used_sites;

    for (auto nn : start_point->next) {
        std::string alt = "";

        while (nn != end_point) {
            if (prg_map.find(nn) != prg_map.end()) {
                auto p_Node = prg_map.at(nn);
                // TODO: The below assertion is a goal we want to achieve.
                //assert(used_sites.find(nn) == used_sites.end());
                used_sites.insert(nn);
                alt += (p_Node->sequence);
                nn = p_Node->next;
            } else {
                try{ // Has this char been committed to PRG string already?
                    if (fixed_point_numbers.at(nn) == 0 ) {;} // Yes: Do nothing
                }
                catch(const std::out_of_range &e) {alt += (nn->characters);} //No: commit the char
                nn = *(nn->next.begin());
            }
        }
        if (alt.length() == 0)  direct_deletion = true;
        alts.push_back(alt);
    }

    num_bubbles_to_process -= 1;
    // Make a prg sequence.
    std::string prg_Seq;

    // Make sure we sort the alleles, so that unit testing is simplified.
    std::sort(alts.begin(), alts.end());
    for (auto s : alts){
        // Prepend the common string to each allele
        if (direct_deletion) s = (start_point->characters) + s;
       prg_Seq+= s + ",";
    }
    prg_Seq.pop_back();
    prg_Seq = "[" + prg_Seq + "]";

    // Prepend the common string pre bifurcation.
    if (start_point->characters != SOURCE_CHAR && !direct_deletion) prg_Seq = start_point->characters + prg_Seq;

    // Postpend the common characters post joining
    // It is OK to postpend at the last incident bubble, because the bubbles are ordered topologically (innermost first)
    if (num_bubbles_to_process == 0 && bubble_map.find(end_point) == bubble_map.end() && end_point->characters != SINK_CHAR){
      prg_Seq += end_point->characters;
    }

    BOOST_LOG_TRIVIAL(debug) << prg_Seq ;

    // Make a prg node, and make it available.
    prg_Node* new_Node = new prg_Node(prg_Seq, end_point);
    prg_map.insert(std::make_pair(start_point,new_Node));
}

void nested_prg::add_site_numbers_and_binary_encode() {
   std::stack<int> marker_stack;
   int max_var_marker{3};
   int char_count{0};

   sdsl::int_vector<> encoded_prg(this->prg.length(), 0, 32);
   BOOST_LOG_TRIVIAL(debug) << "Prg pre serialisation: " << prg;
   for (int i = 0; i<prg.size(); ++i){
       const auto &c = prg[i];

      switch(c) {
          case '[' : {
              max_var_marker += 2;
              marker_stack.push(max_var_marker);
              encoded_prg[char_count++] = max_var_marker;
              break;
          }

          case ']' : {
              assert(!marker_stack.empty());
              encoded_prg[char_count++] = marker_stack.top() + 1;
              marker_stack.pop();
              break;
          }

          case ',' : {
              assert(!marker_stack.empty());
              encoded_prg[char_count++] = marker_stack.top() + 1;
              break;
          }

          default : {
              try {
                  encoded_prg[char_count++] = encode_char(c);
                  break;
              }
              catch(std::exception& e){
                  BOOST_LOG_TRIVIAL(error) << e.what();
                  exit(1);
              }
          }
      }
   }

   this->encoded_prg = encoded_prg;
   BOOST_LOG_TRIVIAL(info) << "Number of sites produced: " << (max_var_marker -3 ) / 2;
}

