#include "FA.hpp"
#include <string>

class prg_Node{

public:
    prg_Node();
    prg_Node(std::string sequence, prg_Node* branch_point);

    friend class nested_prg;

private:
    std::string sequence;
    prg_Node* branch_point;
    bool multifurc;
};

class nested_prg_greedy{

public:
    nested_prg_greedy(auto_Node* root);
    std::string prg;

private:
};
