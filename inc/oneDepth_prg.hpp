#include "FA.hpp"
#include <string>
#include <vector>

class oneDepth_prg {
public:
    std::string prg;
    int num_var_sites;

    oneDepth_prg(auto_Node* root);

private:
    auto_Node* find_fixed_point(auto_Node* starting_point);

    std::vector<std::string> var_region;
    void DFFS(auto_Node* cur_Node, std::string alt, auto_Node* fixed_point);

    std::string serialise_var_region();
};
