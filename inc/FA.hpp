#include "MSA.hpp"
#include <set>
#include <unordered_map>

#ifndef HPP_FA
#define HPP_FA

#define SOURCE_CHAR '#'
#define SINK_CHAR '$'

class auto_Node {
public:
    auto_Node();

    auto_Node(char l, int pos);

    void mark_as_fixed_point() { this->fixed_Point = true; }

    friend class FA;
    friend class oneDepth_prg;
    friend class nested_prg;


private:
    char letter;
    std::set<auto_Node*> next; // Outgoing edges
    std::set<auto_Node*> prev; // Incoming edges
    bool fixed_Point;
    int pos; // Alignment column number.

};

class FA{
public:
    FA(MSA &msa);
    auto_Node* root;
};

#endif
