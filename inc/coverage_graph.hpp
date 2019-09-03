#include "sequence_graph.hpp"

class coverage_Node{
public:
    coverage_Node(){;};
    coverage_Node(const std::string seq, int const pos, int const site_ID = 0, int const allele_ID = 0);
    coverage_Node(auto_Node const& node_in, int const site_ID, int const allele_ID);
    bool in_site() {return is_in_site;};

    friend class coverage_Graph;
private:
    std::string sequence;
    int site_ID;
    int allele_ID;
    int pos;
    std::vector<uint64_t> coverage;
    bool is_in_site;
    std::set<std::shared_ptr<coverage_Node>> prev;
    std::set<std::shared_ptr<coverage_Node>> next;
};

/**
* This class allows to create a copy of a sequence graph. This modified copy is ready to record coverage
* in gramtools quasimap.
* The graph is built 'bubble-up', with children bubbles first. In each bubble,
* each path thus becomes a linear traversal from bubble start to bubble end.
* The nodes are compressed whenever possible, meaning connected nodes both with in and outdegree
* of 1 get joined.
**/
class coverage_Graph{
public:
    /**
     * Build a coverage graph from an existing sequence graph.
     */
    coverage_Graph(sequence_Graph const& graph_in);
private:
    /** Maps the start of a local bubble, to its end.
     * Children nodes appear before parent nodes.
     */
    std::map<std::shared_ptr<coverage_Node>,std::shared_ptr<coverage_Node>,
            std::greater<std::shared_ptr<coverage_Node>> > bubble_map;
};