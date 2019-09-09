#include "MSA.hpp"
#include "utils.hpp"
#include <set>
#include <unordered_map>
#include <memory>

#ifndef HPP_FA
#define HPP_FA

#define SOURCE_CHAR "#"
#define SINK_CHAR "$"

class auto_Node {
    using seqG_ptr = std::shared_ptr<auto_Node>;
public:
    auto_Node();
    ~auto_Node(){
        //std::cout << "Destructor called on " << pos << ":" << characters << std::endl;
        };

    auto_Node(std::string l, int pos);

    void mark_as_fixed_point() { this->fixed_Point = true; }

    /**
     * Compare pointers to `auto_Node`; used in topological ordering (lastmost sequence position first)
     * Equivalence in a set is defined using this, so we also test whether the pointers are the same objects.
     */
    friend bool operator > (const seqG_ptr& lhs, const seqG_ptr& rhs);

    friend std::ostream& operator <<(std::ostream& stream, const seqG_ptr& node);

    friend class FA;
    friend class oneDepth_prg;
    friend class sequence_Graph;
    friend class coverage_Node;
    friend class cov_graph_Constructor;
    template<typename graph_T, typename node_T> friend class stringified_PRG;


private:
    std::string sequence;
    std::set<seqG_ptr> next; // Outgoing edges
    std::set<seqG_ptr> prev; // Incoming edges
    bool fixed_Point;
    int pos; // Alignment column number.

};

using seqG_ptr = std::shared_ptr<auto_Node>;

class FA{
public:
    /**
     * Build Finite Automaton from full multiple-sequence alignment
     */
    FA(MSA& msa);

    /**
     * Build it from delimited boundaries
     * @param haplotypic_resolution the target sequence size of each node.
     */
    FA(MSA& msa, seqG_ptr start_point,
            seqG_ptr end_point, int haplotypic_resolution);
    seqG_ptr root;
    seqG_ptr getSink() const {return sink;};
private:
    seqG_ptr sink;
    const std::set<char> gapping_chars = {'-', '.'};
};

#endif
