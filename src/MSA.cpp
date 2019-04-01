#include <iostream>
#include "MSA.hpp"


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
