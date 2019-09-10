#include "coverage_graph.hpp"

coverage_Node::coverage_Node(std::string const seq, int const pos, int const site_ID, int const allele_ID) :
        sequence(seq), pos(pos), site_ID(site_ID), allele_ID(allele_ID) {
    if (site_ID == 0 && allele_ID == 0) is_in_site = false;
    else is_in_site = true;

    coverage = std::vector<uint64_t>{sequence.length(), 0};
}

coverage_Node::coverage_Node(auto_Node const &node_in, int const site_ID, int const allele_ID) :
        sequence(node_in.sequence), pos(node_in.pos), site_ID(site_ID), allele_ID(allele_ID) {
    if (site_ID == 0 && allele_ID == 0) is_in_site = false;
    else is_in_site = true;

    coverage = std::vector<uint64_t>{sequence.length(), 0};
}

bool operator>(const covG_ptr &lhs, const covG_ptr &rhs) {
    if (lhs->pos == rhs->pos) {
        return lhs.get() > rhs.get();
    } else return lhs->pos > rhs->pos;
}

coverage_Graph::coverage_Graph(sequence_Graph const &graph_in) {
    // Throughout, we must distinguish between nodes in the constructor parameter's graph, and nodes in the new graph.
    cov_graph_Constructor constructor(this, graph_in);

    // Copy & wire all the bubbles
    constructor.parse_bubbles();

    // Copy & wire the root node
    this->root = constructor.make_root(graph_in.root);

    // Copy & wire the non-variant sequence
    constructor.linearTraversal();
}

void cov_graph_Constructor::parse_bubbles() {
    // We will copy each bubble in dependency order: child bubbles first (that is how they are ordered in the map)
    for (auto &s : *in_g_bubble_map) {
        /**
         * Set-up data structures
         */
        allele_ID = 0; // Reset allele ID
        site_ID++;

        // Make entry and exit points;
        // The bubble exit is not initialised here, as it can be known from a previous bubble.
        auto bubble_entry = std::make_shared<coverage_Node>(*s.first, site_ID, allele_ID);
        covG_ptr bubble_exit;

        cur_pos = s.first->pos + 1;
        skip_fixed_point = false;
        deletion_bubble = false;
        deletion_prefix = "";

        // Determine whether we have a direct deletion; if so we need to add some sequence
        // Meaning, we will take the last character of the bubble entry's sequence, take it out from the bubble entry,
        // And preprend it to each allele in the bubble.
        // TODO: do not do this if the bubble starts at SOURCE_CHAR
        for (auto &nn : s.first->next) {
            if (nn == s.second) {
                if (!deletion_bubble) {
                    deletion_bubble = true;
                    deletion_prefix = bubble_entry->sequence.back();
                    bubble_entry->sequence = bubble_entry->sequence.substr(0, bubble_entry->sequence.size() - 1);
                    cur_pos--;
                }
            }
        }

        // Process each allele of the bubble
        bool first_allele{true};
        for (auto seed : s.first->next){
            cur_Node = seed;
            if (first_allele){ // Specify the bubble exit node, once.
                bubble_exit = findBubbleEnd(s.second);
                first_allele = false;
            }
            backWire = bubble_entry;
            parse_allele(s.second);
        }

        // Register the copied bubble
        entry_translation_map.insert(std::make_pair(s.first, bubble_entry));
        bubble_exit = exit_translation_map.at(s.second);
        to_construct->bubble_map.insert(std::make_pair(bubble_entry, bubble_exit));
    }
}

covG_ptr cov_graph_Constructor::parse_allele(seqG_ptr const &bubble_end){
    covG_ptr bubble_exit;
    seqBuffer = "";
    allele_ID++;

    // Check for direct deletion and if present, make nodes with non-empty sequence
    if (deletion_bubble) {
        seqBuffer += deletion_prefix;
        if (cur_Node == bubble_end) {
            auto new_Node = std::make_shared<coverage_Node>(seqBuffer, cur_pos, site_ID, allele_ID);
            wire(new_Node);
            backWire = new_Node;
        }
    }

    // Linear traversal to the end of the allele
    while (cur_Node != bubble_end) {
        if (cur_Node->next.size() == 1) {
            simpleAdvance();
            // If we have reached the bubble end, we need to look at the sequence buffer.
            if (cur_Node != bubble_end) continue;
        }
        if (seqBuffer.size() > 0) make_sequence();
        if (cur_Node == bubble_end) break; //The bubble end could be a bubble_start, and then we would skip it below.
        if (in_g_bubble_map->find(cur_Node) != in_g_bubble_map->end()) bubble_skip();
    }

    bubble_exit = exit_translation_map.at(bubble_end);

    // The backWire can be bubble_exit if the allele path ends by going through a bubble ending at bubble exit.
    // If that's the case, we want to avoid self linking.
    if (backWire != bubble_exit) wire(bubble_exit);
}


void cov_graph_Constructor::linearTraversal() {
    site_ID, allele_ID = 0, 0; // Any sequence made in this function will be outside variant sites.
    seqBuffer = "";

    while (cur_Node->sequence != SINK_CHAR) {
        if (cur_Node->next.size() == 1) {
            simpleAdvance();
            if (cur_Node->sequence != SINK_CHAR) continue;
        }
        if (seqBuffer.size() > 0) make_sequence();
        // Case: we have a bubble
        if (cur_Node->sequence != SINK_CHAR) bubble_skip();
    }

    assert(cur_Node->sequence == SINK_CHAR);
    if (backWire->sequence != SINK_CHAR) { // The sink is not a bubble end, it does not exist yet.
        auto new_Node = std::make_shared<coverage_Node>(SINK_CHAR, cur_Node->pos, 0, 0);
        wire(new_Node);
    }
}


void cov_graph_Constructor::make_sequence() {
    auto new_Node = std::make_shared<coverage_Node>(seqBuffer, cur_pos, site_ID, allele_ID);
    wire(new_Node);
    backWire = new_Node;
    seqBuffer = "";
}

void cov_graph_Constructor::bubble_skip(){
    auto &translated_bubble = entry_translation_map.at(cur_Node); // Will throw error if not there; it ought to be.
    // The backWire will be the same as the bubble entry if there is a
    // node which is both a bubble end and a bubble start. In that case, no link needed.
    if (backWire != translated_bubble) wire(translated_bubble);

    backWire = to_construct->bubble_map.at(translated_bubble);
    cur_Node = in_g_bubble_map->at(cur_Node);
    cur_pos = cur_Node->pos + 1;
}

void cov_graph_Constructor::simpleAdvance() {
    // Only commit the sequence if we are not in a bubble end.
    if (exit_translation_map.find(cur_Node) == exit_translation_map.end()) {
        seqBuffer += cur_Node->sequence;
    }
    cur_Node = *(cur_Node->next.begin());
}


covG_ptr cov_graph_Constructor::make_root(seqG_ptr const &in_root){
    covG_ptr root;
    // copy the root node sequence
    cur_Node = in_root;
    assert(cur_Node->sequence == SOURCE_CHAR);
    if (cur_Node->next.size() > 1) { // Case: the root is a bubble start
        auto &translated_bubble = entry_translation_map.at(
                cur_Node); // Will throw error if not there; it ought to be.
        backWire = to_construct->bubble_map.at(translated_bubble);
        cur_Node = in_g_bubble_map->at(cur_Node);
        cur_pos = cur_Node->pos + 1;
        root = translated_bubble;
    } else {
        cur_pos = cur_Node->pos;
        auto new_Node = std::make_shared<coverage_Node>(SOURCE_CHAR, cur_pos, 0, 0);
        backWire = new_Node;
        cur_Node = *(cur_Node->next.begin());
        root = new_Node;
    }
    return root;
}

void cov_graph_Constructor::wire(covG_ptr& target){
    backWire->next.push_back(target);
    target->prev.push_back(backWire);
}

covG_ptr cov_graph_Constructor::findBubbleEnd(seqG_ptr const &bubble_end){
    covG_ptr bubble_exit;
    // We might end in a bubble that we have ended in before.
    // If we do, reuse it!
    if (exit_translation_map.find(bubble_end) != exit_translation_map.end())
        bubble_exit = exit_translation_map.at(bubble_end);
    else { // Never seen it before: make it
        if (entry_translation_map.find(bubble_end)
        != entry_translation_map.end()) { // Case: the bubble ends where a bubble starts; use that
            auto translated_bubble = entry_translation_map.at(bubble_end);
            bubble_exit = translated_bubble;
        } else bubble_exit = std::make_shared<coverage_Node>(*bubble_end, 0, 0);

        exit_translation_map.insert({bubble_end, bubble_exit});
        // Map the number of incidents to the bubble if first time we look at bubble end.
        int num_incidents = fixed_point_numbers->at(bubble_end);
        to_construct->fixed_point_numbers.insert(std::make_pair(bubble_exit, num_incidents));
    }
    return bubble_exit;
}
