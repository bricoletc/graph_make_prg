#include "FA.hpp"
#include <string>
#include <vector>

class oneDepth_prg {
public:
    std::string prg;
    int num_var_sites;

    oneDepth_prg(std::shared_ptr<auto_Node> root);

private:
    std::shared_ptr<auto_Node> find_fixed_point(std::shared_ptr<auto_Node> starting_point);

    std::vector<std::string> var_region;
    void DFFS(std::shared_ptr<auto_Node> cur_Node, std::string alt, std::shared_ptr<auto_Node> fixed_point);

    std::string serialise_var_region();
};
