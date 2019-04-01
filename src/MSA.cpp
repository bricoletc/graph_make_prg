#include <iostream>
#include "MSA.hpp"


MSA::MSA(std::string MSA_fname) {
    num_records = 0;
    offset = 0;

    handle = std::ifstream(MSA_fname, std::ifstream::binary);
    if (handle.fail()) {
        std::cout << "Error: cannot open file " << MSA_fname << ".\n Exiting.";
        exit(1);
    }
    find_starts();

    more_columns = true;
}


MSA::~MSA() {
    handle.close();
}


// TODO: add check that each record has same number of chars
void MSA::find_starts() {

    char c = '\0';

    while (!handle.eof()) {

        handle.get(c);
        if (c != '>') continue;

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
