#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#ifndef HPP_MSA
#define HPP_MSA

class MSA {
public:
    MSA(std::string MSA, bool is_file = true);
    ~MSA();
    /**
     * Extracts the next unprocessed character for each Fasta record
     * @return
     */
    std::vector<char> next_column();
    /**
     * @param new_offset the alignment column number from which the next characters will come.
     */
    void reposition(int new_offset);
    int num_records;
    bool more_columns;
    int record_size;


private:
    bool is_a_file;
    std::stringstream ss_handle;
    std::ifstream if_handle;
    /**
     * stringstream and ifstream are both derived classes from istream
     * This function takes either and returns a uniform object on which to operate regardless of MSA data input
     */
    std::istream& get_handle();
    std::vector<int> seekg_starts;
    int offset; // Stores the current alignment column number (zero-based).
    void find_starts(std::istream& handle);
};

#endif
