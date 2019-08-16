#include "oneDepth_prg.hpp"
#include <unordered_map>

oneDepth_prg::oneDepth_prg(std::shared_ptr<auto_Node> root) {
    std::shared_ptr<auto_Node> cur_Node = root;
    prg = "";
    num_var_sites = 0;


    while (cur_Node->letter != SINK_CHAR) {
        if (cur_Node->next.size() == 1) {
            std::shared_ptr<auto_Node> nn = *(cur_Node->next.begin());

            if (nn->letter != SINK_CHAR) prg += nn->letter;

            cur_Node = nn;
        } else {
            std::shared_ptr<auto_Node> fixed_point = NULL;
            std::string alt = "";

            // Look for direct link to a fixed point
            for (auto node : cur_Node->next) {
                if (node->fixed_Point) {
                    fixed_point = node;
                    alt = prg[prg.size() - 1];
                    prg.pop_back();
                    break;
                }
            }

            if (fixed_point == NULL) fixed_point = find_fixed_point(cur_Node);

            var_region.clear();
            DFFS(cur_Node, alt, fixed_point); //populate var_region
            prg += serialise_var_region();

            num_var_sites++;

            cur_Node = fixed_point;


            if (cur_Node->letter != SINK_CHAR) prg += cur_Node->letter;
        }
    }
};

/**
 * DFS traversal from the given node until we find a fixed point.
 */
std::shared_ptr<auto_Node> oneDepth_prg::find_fixed_point(std::shared_ptr<auto_Node> cur_Node) {
    //If the current node is already a fixed point, avoid returning it.
    if (cur_Node->fixed_Point) cur_Node = *(cur_Node->next.begin());

    while (!cur_Node->fixed_Point) {
        for (auto node : cur_Node->next){
            cur_Node = node;
            break;
        }
    }
    return cur_Node;
}

/**
 * Depth-First Fixed-point bounded Search.
 */
void oneDepth_prg::DFFS(std::shared_ptr<auto_Node> cur_Node, std::string alt, std::shared_ptr<auto_Node> fixed_point) {
    for (auto node : cur_Node->next) {
        if (node == fixed_point){
            var_region.push_back(alt);
        }
        else {
            std::string new_alt = alt + node->letter;
            DFFS(node, new_alt, fixed_point);
        }
    }
}

std::string oneDepth_prg::serialise_var_region() {
    std::string serialised = "";

    int var_site_num = num_var_sites * 2 + 5;
    serialised += std::to_string(var_site_num);

    std::string allele;

    for (int i = 0; i < var_region.size(); i++) {
        allele = var_region[i];

        // Add allele marker
        if (i > 0) serialised += std::to_string(var_site_num + 1);

        // Add allele
        serialised += allele;
    }
    serialised += std::to_string(var_site_num);

    return serialised;
}
