#include "FA.hpp"
#include "utils.hpp"
#include <unordered_map>
#include <stack>
#include <queue>

#ifndef HPP_NESTED_PRG
#define HPP_NESTED_PRG


struct incidence_fixed_point{
    seqG_ptr fixed_point;
    int num_incidents;
    int haplotype_resolution;
    seqG_ptr earliest_incident;
    int pos_earliest_incident;

    /**
     * For ordering STLs containing this object
     */
    friend bool operator < (const incidence_fixed_point& lhs, const incidence_fixed_point& rhs);
};


class sequence_Graph{
public:
    sequence_Graph(seqG_ptr root, std::string MSA_file = "", bool is_file = true, int max_num_incidents = 2);

    friend class coverage_Graph;
    friend class cov_graph_Constructor;
    template<typename graph_T, typename node_T> friend class stringified_PRG;

private:
    /**
     * Data structures
     */
     seqG_ptr root;

    // User definable parameter controlling how much recombination we want to allow between haplotypes.
    int max_num_incidents;

    // Path to file containing MSA that can be loaded to rewrite portions of graph.
    std::string MSA_file;
    bool is_file;

    /** Maps the start of a local bubble, to its end.
     * Children nodes appear before parent nodes.
     * */
    std::map<seqG_ptr,seqG_ptr,
        std::greater<seqG_ptr> > bubble_map;



    /** Maps the end of a local bubble to the number of bubbles ending in it left to process. When that number is 0,
     * the corresponding character can be committed to the prg string.*/
    std::unordered_map<seqG_ptr,int> fixed_point_numbers;

    /**
     * Maps the end of a local bubble to a struct containing:
     *  - The earliest (positionally) bubble that end into it.
     *  - The number of bubbles that end into it.
     */
    std::unordered_map<seqG_ptr,incidence_fixed_point> fixed_point_incidence_map;

    std::set<incidence_fixed_point,
        std::less<incidence_fixed_point>> large_incidence_fixed_points;


    /**
     * Functions
     */

    /**
     * Traverses full NFA graph and maps all bubbles.
     * @see map_bubbles()
     */
    void map_all_bubbles(seqG_ptr root);

    /**
     * Finds the first `auto_Node` at which all paths going from `start_point` end.
     */
    seqG_ptr map_bubbles(seqG_ptr start_point, int haplotype_res);

    /**
     * Rewrite portions of the graph with decreasing levels of recombination.
     */
    void haplotype_expand_bubbles();

    /**
     * Free the heap allocated nodes in between (and not including) the two specified pointers
     * Specify the graph source and sink nodes to free the whole graph.
     */
    void delete_in_between(seqG_ptr start_point, seqG_ptr end_point);

    int rebuild_in_between(MSA& msa, incidence_fixed_point& i);

    /**
     * Populate a set with high incidence fixed points
     */
    void populate_large_incidences();


};

#endif
