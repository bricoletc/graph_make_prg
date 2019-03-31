#include <iostream>
#include <unordered_map>
#include "msa_to_dfa.hpp"

auto_Node::auto_Node()
        :
        letter('\0'),
        fixed_Point(false) {
}

auto_Node::auto_Node(char l)
        :
        letter(l),
        fixed_Point(false) {
}


MSA::MSA(std::string MSA_fname) {
    num_records = 0;
    offset = 0;

    fhandle = std::ifstream(MSA_fname, std::ifstream::binary);
    if (fhandle.fail()) {
        std::cout << "Error: cannot open file " << MSA_fname << ".\n Exiting.";
        exit(1);
    }
    find_starts();

    more_columns = true;
}


MSA::~MSA() {
    fhandle.close();
}


// TODO: add check that each record has same number of chars
void MSA::find_starts() {

    char c = '\0';

    while (!fhandle.eof()) {

        fhandle.get(c);
        if (c != '>') continue;

        // Go to end of header line
        while (c != '\n' && !fhandle.eof()) fhandle.get(c);


        // Make sure there is at least one valid char
        char t = fhandle.peek();
        if (t == '\n') {
            std::cout << "Found a record with newline between fasta header and sequence. Record will be ignored"
                      << std::endl;
            continue;
        }

        if (fhandle.eof()) {
            std::cout << "Found a header only fasta record. Record will be ignored." << std::endl;
            continue;
        }

        // Seekg needs to give position before a .get() for the first column element.
        seekg_starts.push_back(fhandle.tellg());

        num_records++;
    }


    if (c == '\0') {
        std::cout << "Error: the MSA file provided does not contain any fasta header" << std::endl;
        std::exit(1);
    }

    if (seekg_starts.size() == 0) {
        std::cout << "Error: cannot find any valid records in MSA file provided" << std::endl;
        std::exit(1);
    }
}


std::vector<char> MSA::next_column() {
    std::vector<char> column;

    char c;
    fhandle.clear();

    for (int i = 0; i < num_records; i++) {
        int s = seekg_starts[i];

        fhandle.seekg(s + offset, std::ios::beg);

        fhandle.get(c);

        if (c == '\n') fhandle.get(c);

        if (c == '\n') {
            std::printf("Error: empty line found inside record at byte %i. Exiting.", s + offset);
            exit(1);
        }

        column.push_back(c);

        // Specific to the last record. Look for eof and signal there are no more columns if found.
        if (i == num_records - 1) {
            char t;
            fhandle.get(t);
            while (t == '\n' && !fhandle.eof()) fhandle.get(t);

            if (fhandle.eof()) more_columns = false;
        }

    }

    offset++;

    return column;
}


FA::FA(MSA &msa) {
    static std::set<char> gapping_chars = {'-', '.'};

    root = auto_Node(SOURCE_CHAR);
    int N = msa.num_records;

    auto_Node *cur_Nodes[N];

    // Initialise to point to root.
    for (int i = 0; i < N; i++) {
        cur_Nodes[i] = &root;
    }

    std::vector<char> letters;

    std::unordered_map<char, int> new_Nodes;

    while (msa.more_columns) {
        letters = msa.next_column();
        for (int i = 0; i < N; i++) {
            char l = letters[i];

            // If character indicating indel padding is found, do nothing for this letter.
            if (gapping_chars.find(l) != gapping_chars.end()) continue;

            // Create a new node; and make it accessible.
            if (new_Nodes.find(l) == new_Nodes.end()){
                auto_Node new_Node = auto_Node(l);
                all_Nodes.push_back(new_Node);
                new_Nodes.insert(std::make_pair(l, all_Nodes.size() - 1));
            }

            auto_Node ret_node = all_Nodes[new_Nodes.at(l)];
            auto_Node prec_node = *(cur_Nodes[i]);

            // If there is already an edge between these two nodes, carry on.
            if (prec_node.has_edgeTo(&ret_node)) continue;

            prec_node.make_edgeTo(&ret_node);

            // Make retrieved node accessible to next column.
            cur_Nodes[i] = &ret_node;

        }

        if (new_Nodes.size() == 1) cur_Nodes[0]->mark_as_fixed_point();

        new_Nodes.clear();
    }

    // Final link to a fixed_point auto_Node
    auto_Node sink_node = auto_Node(SINK_CHAR);
    all_Nodes.push_back(sink_node);

    for (int i = 0; i < N; i++) {
        auto_Node cn = *(cur_Nodes[i]);

        if (not cn.has_edgeTo(&sink_node)) cn.make_edgeTo(&all_Nodes[all_Nodes.size() - 1]);
    }
    sink_node.mark_as_fixed_point();

}


oneDepth_prg::oneDepth_prg(auto_Node &root) {
    auto_Node *cur_Node = &root;
    prg = "";
    num_var_sites = 0;


    while (cur_Node->letter != SINK_CHAR) {
        if (cur_Node->next.size() == 1) {
            auto_Node *nn = *(cur_Node->next.begin());

            if (nn->letter != SINK_CHAR) prg += nn->letter;

            cur_Node = nn;
        } else {
            auto_Node *fixed_point = NULL;
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
auto_Node *oneDepth_prg::find_fixed_point(auto_Node *cur_Node) {
    // Test starting condition; the source node is a fixed point but we do not want to return it.
    if (cur_Node->letter == SOURCE_CHAR) cur_Node = *(cur_Node->next.begin());

    while (!cur_Node->fixed_Point) {
        cur_Node = *(cur_Node->next.begin());
    }
    return cur_Node;
}

std::vector<std::string> oneDepth_prg::DFFS(auto_Node *cur_Node, std::string alt, auto_Node *fixed_point) {
    for (auto node : cur_Node->next) {
        if (node == fixed_point) var_region.push_back(alt);
        else DFFS(node, alt + node->letter, fixed_point);
    }
}

std::string oneDepth_prg::serialise_var_region() {
    std::string serialised = "";

    int var_site_num = num_var_sites * 2 + 4;
    serialised += std::to_string(var_site_num);

    std::string allele;

    for (int i = 0; i < var_region.size(); i++) {
        allele = var_region[i];

        // Add allele marker
        if (i > 0) serialised += std::to_string(var_site_num + 1);

        // Add allele
        serialised += allele;
    }

    return serialised;
}
