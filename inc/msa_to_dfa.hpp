#include <set>
#include <string>
#include <vector>
#include <fstream>

class auto_Node {
public:
	auto_Node();

	auto_Node(char l);


	bool has_edgeTo(auto_Node n) {
		std::set<auto_Node*> edges = this->next;

		return edges.find(&n) != edges.end();
	}

	void make_edgeTo(auto_Node n) { this->next.insert(&n); }

	void mark_as_fixed_point() { this->fixed_Point = true; }


private:
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


auto_Node MSA_to_FA(MSA &msa);
