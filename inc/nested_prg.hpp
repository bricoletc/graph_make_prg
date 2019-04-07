#include "FA.hpp"
#include <unordered_map>
#include <stack>

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

private:
    std::unordered_map<auto_Node*,auto_Node*> bubble_map;
    std::stack<auto_Node*> topological_order;

    std::unordered_map<auto_Node*, prg_Node*> prg_map;
    std::unordered_map<auto_Node*,std::set<auto_Node*>> fixed_point_map;

    void map_all_bubbles(auto_Node* root);
    auto_Node* map_bubbles(auto_Node* start_point);
    void parse_bubbles(auto_Node* start_point, auto_Node* end_point);

};
