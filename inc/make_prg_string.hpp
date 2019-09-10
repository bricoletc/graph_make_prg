#include "sequence_graph.hpp"
#include "coverage_graph.hpp"

template<typename T>
class prg_Node{

public:
    prg_Node();
    prg_Node(std::vector<uint64_t> sequence, std::shared_ptr<T> next);

    std::vector<uint64_t> sequence;
    std::shared_ptr<T> next;
};

/**
 * This class builds a PRG string from a directed, acyclic graph of sequence.
 * It can do this from a sequence graph and a coverage graph both.
 */
template<typename graph_T, typename node_T>
class stringified_PRG{
public:
    stringified_PRG(graph_T& graph_in);
    std::vector<uint64_t> prg;

private:
    /** Allows checking if a sequence node has a `prg_Node`, in which case the `prg_Node` will be used.*/
    std::unordered_map<std::shared_ptr<node_T>, prg_Node<node_T>*> prg_map;

    // A local reference to the in graph's structure
    std::map<std::shared_ptr<node_T>, std::shared_ptr<node_T>,
    std::greater<std::shared_ptr<node_T> > > bubble_map;

    // A local reference to the in graph's structure
    std::unordered_map<std::shared_ptr<node_T>, int> fixed_point_numbers;
    /**
     * Traverses a bubble using a depth first search; builds all alleles and wraps them in a prg string.
     */
    void parse_bubble(std::shared_ptr<node_T> start_point, std::shared_ptr<node_T> end_point, int site_ID);
};


std::vector<uint64_t> prg_string_to_ints(std::string const& string_prg);
