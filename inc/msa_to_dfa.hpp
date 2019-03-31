#include <set>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#define SOURCE_CHAR '#'
#define SINK_CHAR '$'

class auto_Node {
public:
	auto_Node();

	auto_Node(char l);


	bool has_edgeTo(auto_Node* n) {
		std::set<auto_Node*> edges = this->next;

		return edges.find(n) != edges.end();
	}

	void make_edgeTo(auto_Node* n) {
		this->next.insert(n);
	}

	void mark_as_fixed_point() { this->fixed_Point = true; }

	friend class oneDepth_prg;


//private:
	char letter;
	std::set<auto_Node*> next;
	bool fixed_Point;

};

class MSA {
public:
	MSA(std::string MSA_fname);
	~MSA();
	std::vector<char> next_column();
	int num_records;
	bool more_columns;


private:
	std::ifstream fhandle;
	std::vector<int> seekg_starts;
	int offset;
	void find_starts();
};


class FA{
public:
	FA(MSA &msa);
	auto_Node root;

//private:
	std::vector<auto_Node> all_Nodes;
};


class oneDepth_prg {
public:
	std::string prg;
	int num_var_sites;

	oneDepth_prg(auto_Node &root);

private:
	auto_Node* find_fixed_point(auto_Node* starting_point);

	std::vector<std::string> var_region;
	std::vector<std::string> DFFS(auto_Node* cur_Node, std::string alt, auto_Node* fixed_point);

	std::string serialise_var_region();
};
