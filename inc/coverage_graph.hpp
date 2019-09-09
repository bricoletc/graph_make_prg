#include "sequence_graph.hpp"


class coverage_Node{
    using covG_ptr = std::shared_ptr<coverage_Node>;
public:
    coverage_Node(){;};
    coverage_Node(const std::string seq, int const pos, int const site_ID = 0, int const allele_ID = 0);
    coverage_Node(auto_Node const& node_in, int const site_ID, int const allele_ID);
    bool in_site() {return is_in_site;};

    /**
     * Compare pointers to `coverage_Node`; used in topological ordering (lastmost sequence position first)
     * Equivalence in a set is defined using this, so we also test whether the pointers are the same objects.
     */
    friend bool operator > (const covG_ptr& lhs, const covG_ptr& rhs);

    friend class coverage_Graph;
    friend class cov_graph_Constructor;
    template<typename graph_T, typename node_T> friend class stringified_PRG;

private:
    std::string sequence;
    int site_ID;
    int allele_ID;
    int pos;
    std::vector<uint64_t> coverage;
    bool is_in_site;
    std::set<covG_ptr> prev;
    std::set<covG_ptr> next;
};

using covG_ptr = std::shared_ptr<coverage_Node>;
using seqG_map = const std::map<seqG_ptr, seqG_ptr, std::greater<seqG_ptr>>;

/**
* This class allows to create a copy of a sequence graph. This modified copy is ready to record coverage
* in gramtools quasimap.
* The graph is built 'bubble-up', with children bubbles first. In each bubble,
* each path thus becomes a linear traversal from bubble start to bubble end.
**/
class coverage_Graph{
public:

    covG_ptr root;

    /**
     * Build a coverage graph from an existing sequence graph.
     */
    coverage_Graph(sequence_Graph const& graph_in);

    /** Maps the start of a local bubble, to its end.
     * Children nodes appear before parent nodes.
     */
    std::map<covG_ptr,covG_ptr, std::greater<covG_ptr> > bubble_map;

    /** Maps a bubble end to the number of bubbles that end there.
     */
     std::unordered_map<covG_ptr, int> fixed_point_numbers;

     friend class cov_graph_Constructor;
};

/**
 * A class in charge of copying a sequence graph into a coverage graph
 *
 * Features:
 * - The nodes are compressed whenever possible, meaning connected nodes both with in and outdegree of 1 get joined.
 * - Bubbles where there is a direct edge between bubble start and bubble end nodes get 'expanded': each allele gets
 * the bubble start's last character as a prefix. This later avoids empty allelic sequences.
 */
class cov_graph_Constructor{
public:
    cov_graph_Constructor(coverage_Graph* to_construct, sequence_Graph const& graph_in) :
        site_ID{0},
        to_construct{to_construct},
        in_g_bubble_map{&(graph_in.bubble_map)},
        fixed_point_numbers{&(graph_in.fixed_point_numbers)}{;};

    coverage_Graph* to_construct; // The target graph which will be built.

    /**
     * Functions used to build the cov graph
     */
    void parse_bubbles();

    covG_ptr parse_allele(seqG_ptr const& bubble_end);

    /**
     * Traverses the full sequence graph, copying & wiring sequence outside variant sites.
     */
    void linearTraversal();

    /**
     * Build a new sequence node
     */
    void make_sequence();
    /**
     * Skip a previously copied bubble, with side effects.
     */
    void bubble_skip();
    /**
     * Move directly from node with outdegree one to node with indegree one, with side effects.
     */
    void simpleAdvance();

    /**
     * Build & wire the coverage graph's root node
     */
    covG_ptr make_root(seqG_ptr const &in_root);

    /**
     * Build edge between the most recent node and the @param target
     */
    void wire(covG_ptr& target);

    /**
     * Determine bubble end, either by creating it or re-using it
     */
    covG_ptr findBubbleEnd(seqG_ptr const &bubble_end);

private:
    /**
     * Globally used data
     */
    // maps to the coverage graph bubble start & end equivalents.
    std::unordered_map<seqG_ptr, covG_ptr > entry_translation_map;
    std::unordered_map<seqG_ptr, covG_ptr > exit_translation_map;

    covG_ptr backWire; // Points to latest previous node
    int cur_pos; std::string seqBuffer{""}; // For giving new nodes sequence

    seqG_ptr cur_Node; // For traversing the input graph.
    int site_ID;
    int allele_ID;

    /**
     * Bubble-processing specific data
     **/
    bool skip_fixed_point;
    bool deletion_bubble;
    std::string deletion_prefix;

    /**
     * Input graph specific data
     */
     seqG_map* in_g_bubble_map;

     const std::unordered_map<seqG_ptr, int>* fixed_point_numbers;
};
