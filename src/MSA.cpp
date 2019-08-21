#include <iostream>
#include <sstream>
#include "assert.h"
#include "MSA.hpp"
#include "utils.hpp"

MSA::MSA(std::string MSA, bool is_file) {
    num_records = 0;
    offset = 0;
    is_a_file = is_file;

    if (is_file){
        if_handle = std::ifstream(MSA, std::ifstream::binary);
        if (if_handle.fail()){
            std::cout << "Error: cannot open file " << MSA << ".\n Exiting.";
            exit(1);
        };
        find_starts(if_handle);
    }
    else {
        ss_handle = std::stringstream(MSA);
        if (ss_handle.fail()){
            std::cout << "Error: cannot open provided string as stream .\n Exiting.";
            exit(1);
        }
        find_starts(ss_handle);
    }

    more_columns = true;
}


MSA::~MSA() {
    if (is_a_file) if_handle.close();
}


void MSA::find_starts(std::istream& handle) {

    // We will tally fasta record sizes and abort if they are not all the same.
    bool master_record = true;
    bool in_record = false;
    int master_record_size{0};
    int cur_record_size;

    char c = '\0';

    while (!handle.eof()) {

        handle.get(c);
        if (c != '>'){
            if (in_record && c != '\n') ++cur_record_size;
            continue;
        }

        // Go to end of header line
        while (c != '\n' && !handle.eof()) handle.get(c);


        // Make sure there is at least one valid char
        char t = handle.peek();
        if (t == '\n') {
            std::cout << "Found a record with newline between fasta header and sequence. Record will be ignored"
                      << std::endl;
            continue;
        }

        if (handle.eof()) {
            std::cout << "Found a header only fasta record. Record will be ignored." << std::endl;
            continue;
        }

        // Seekg needs to give position before a .get() for the first column element.
        seekg_starts.push_back(handle.tellg());
        num_records++;

        if (!in_record)
            in_record = true;
        else if (master_record){
            master_record_size = cur_record_size;
            master_record = false;
        }
        else {
            if (cur_record_size != master_record_size){
                BOOST_LOG_TRIVIAL(error) << "Found two records with different "
                                            "lengths in input file: \n"
                                            "record 1 " << master_record_size <<
                                            " vs record " << num_records - 1 <<
                                            " " << cur_record_size << std::endl <<
                                            "Please reformat.";
                exit(1);
            }
        }
        cur_record_size = 0;
    }


    if (c == '\0') {
        std::cout << "Error: the MSA file provided does not contain any fasta header" << std::endl;
        std::exit(1);
    }

    if (seekg_starts.size() == 0) {
        std::cout << "Error: cannot find any valid records in MSA file provided" << std::endl;
        std::exit(1);
    }
    record_size = master_record_size;
}

std::istream& MSA::get_handle(){
   if (is_a_file) return if_handle;
   else return ss_handle;
};

std::vector<char> MSA::next_column() {
    std::istream& handle = get_handle();
    std::vector<char> column;

    char c;
    handle.clear();

    for (int i = 0; i < num_records; i++) {
        int s = seekg_starts[i];

        handle.seekg(s + offset, std::ios::beg);

        handle.get(c);

        if (c == '\n') handle.get(c);

        if (c == '\n') {
            std::printf("Error: empty line found inside record at byte %i. Exiting.", s + offset);
            exit(1);
        }

        column.push_back(c);

        // Specific to the last record. Look for eof and signal there are no more columns if found.
        if (i == num_records - 1) {
            char t;
            handle.get(t);
            while (t == '\n' && !handle.eof()) handle.get(t);

            if (handle.eof()) more_columns = false;
        }

    }
    offset++;

    return column;
}

void MSA::reposition(int new_offset){
    offset = new_offset;
}
