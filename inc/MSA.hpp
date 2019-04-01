#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#ifndef HPP_MSA
#define HPP_MSA

class MSA {
public:
    MSA(std::string MSA_fname);
    ~MSA();
    std::vector<char> next_column();
    int num_records;
    bool more_columns;


private:
    std::ifstream handle;
    std::vector<int> seekg_starts;
    int offset;
    void find_starts();
};

#endif
