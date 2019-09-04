#include "sequence_graph.hpp"
#include "coverage_graph.hpp"
#include <sdsl/suffix_arrays.hpp>

template<typename T>
class prg_Node{

public:
    prg_Node();
    prg_Node(std::vector<uint64_t> sequence, std::shared_ptr<auto_Node> next);


private:
    std::vector<uint64_t> sequence;
    std::shared_ptr<T> next;
};

template<typename T>
class stringified_PRG{
public:
    template<typename T>
    stringified_PRG(T& graph_in);
    std::vector<uint64_t> prg;

private:
    /** Allows checking if an `auto_Node` has a `prg_Node`, in which case the `prg_Node` will be used.*/
    std::unordered_map<std::shared_ptr<T>, prg_Node*> prg_map;

    // A local reference to the in graph's structure
    const std::unordered_map<std::shared_ptr<T>, std::shared_ptr<T>> bubble_map;

    // A local reference to the in graph's structure
    const std::unordered_map<std::shared_ptr<T>, int> fixed_point_numbers;
    /**
     * Traverses a bubble using a depth first search; builds all alleles and wraps them in a prg string.
     */
        template<typename T>
        void parse_bubble(std::shared_ptr<T> start_point, std::shared_ptr<T> end_point, int site_ID);
};


void binary_encode();
