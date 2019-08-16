#include "FA.hpp"
#include "utils.hpp"
#include <unordered_map>
#include <stack>
#include <queue>

#ifndef HPP_NESTED_PRG
#define HPP_NESTED_PRG

class prg_Node{

public:
    prg_Node();
    prg_Node(std::string sequence, auto_Node* next);

    friend class nested_prg;

private:
    std::string sequence;
    auto_Node* next;
};

class nested_prg{
public:
    nested_prg(auto_Node* root);
    std::string prg;
    std::string serialised_prg;

private:
    std::unordered_map<auto_Node*,auto_Node*> bubble_map; /**< Maps the start of a local bubble, to its end. */

    /**< Topological ordering of the FA object; children nodes appear before parent nodes. */
    std::priority_queue<auto_Node *, std::vector<auto_Node *>, std::less<std::vector<auto_Node*>::value_type>> topological_order;


    /** Allows checking if an `auto_Node` has a `prg_Node`, in which case the `prg_Node` will be used.*/
    std::unordered_map<auto_Node*, prg_Node*> prg_map;

    /** Maps the end of a local bubble to a set of direct ancestors. When that set is empty, the corresponding
     * character can be committed to the prg string.*/
    std::unordered_map<auto_Node*,std::set<auto_Node*>> fixed_point_map;

    /**
     * Maps the end of a local bubble to a pair of:
     *  - The earliest (positionally) bubble that end into it.
     *  - The number of bubbles that end into it.
     */
    std::unordered_map<auto_Node*,std::pair<auto_Node*,int>> fixed_point_numbers;

    /**
     * Traverses full NFA graph and maps all bubbles.
     * @see map_bubbles()
     */
    void map_all_bubbles(auto_Node* root);

    /**
     * Finds the first `auto_Node` at which all paths going from `start_point` end.
     */
    auto_Node* map_bubbles(auto_Node* start_point);

    /**
     * Traverses a bubble using a depth first search; builds all alleles and wraps them in a prg string.
     */
    void parse_bubbles(auto_Node* start_point, auto_Node* end_point);

    void serialise_prg();
};

#endif
