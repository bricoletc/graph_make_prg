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
    friend bool operator > (const std::shared_ptr<auto_Node>& lhs, const std::shared_ptr<auto_Node>& rhs);

    friend std::ostream& operator <<(std::ostream& stream, const std::shared_ptr<auto_Node>& node);

    friend class FA;
    friend class oneDepth_prg;
    friend class sequence_Graph;
    friend class coverage_Node;
    friend class coverage_Graph;


private:
    std::string characters;
    std::set<std::shared_ptr<auto_Node>> next; // Outgoing edges
    std::set<std::shared_ptr<auto_Node>> prev; // Incoming edges
    bool fixed_Point;
    int pos; // Alignment column number.

};

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
    FA(MSA& msa, std::shared_ptr<auto_Node> start_point,
            std::shared_ptr<auto_Node> end_point, int haplotypic_resolution);
    std::shared_ptr<auto_Node> root;
    std::shared_ptr<auto_Node> getSink() const {return sink;};
private:
    std::shared_ptr<auto_Node> sink;
    const std::set<char> gapping_chars = {'-', '.'};
};

#endif
