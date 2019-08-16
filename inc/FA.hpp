#include "MSA.hpp"
#include "utils.hpp"
#include <set>
#include <unordered_map>
#include <memory>

#ifndef HPP_FA
#define HPP_FA

#define SOURCE_CHAR '#'
#define SINK_CHAR '$'

class auto_Node {
public:
    auto_Node();
    ~auto_Node(){std::cout << "Destructor called on " << pos << ":" << letter << std::endl;};

    auto_Node(char l, int pos);

    void mark_as_fixed_point() { this->fixed_Point = true; }

    friend bool operator > (const auto_Node& lhs, const auto_Node& rhs);

    friend class FA;
    friend class oneDepth_prg;
    friend class nested_prg;


private:
    char letter;
    std::set<std::shared_ptr<auto_Node>> next; // Outgoing edges
    std::set<std::shared_ptr<auto_Node>> prev; // Incoming edges
    bool fixed_Point;
    int pos; // Alignment column number.

};

class FA{
public:
    FA(MSA &msa);
    std::shared_ptr<auto_Node> root;
};

#endif
