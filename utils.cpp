#include "utils.h"

// read sample name and sequence from fasta file
std::tuple<std::string, std::string> next_sample(std::ifstream &in) {
    std::string name, seq;
    std::string line;

    while (std::getline(in, line)) {
        if (line[0] == '>') {
            name = line.substr(1);
            break;
        }
    }
    while (std::getline(in, line)) {
        if (line[0] == '>') {
            in.seekg(-line.size(), std::ios_base::cur);
            break;
        }
        if (line.size() > 1 && line.back() == '\r') line.pop_back(); // remove newline character
        seq += line;
    }
    // step back one char if we are not at the end of the file
    if (!in.eof()) {
       in.seekg(-1, std::ios_base::cur);
    }

    return std::make_tuple(name, seq);
}

int is_unknown_site(char base, std::string seq_type) {
    if (seq_type == NT) {
        switch (base) {
            case '-':
            case '?':
            case '*':
            case 'X':
            case 'x':
            case 'N':
            case 'n':
                return 1;
            default:
                return 0;
        }
    } else {
        switch (base) {
            case '-':
            case '?':
            case '*':
            case 'X':
            case 'x':
                return 1;
            default:
                return 0;
        }
    }
}
