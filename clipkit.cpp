#include <string>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>
#include <map>

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
        seq += line;
    }
    // step back one char if we are not at the end of the file
    if (!in.eof()) {
       in.seekg(-1, std::ios_base::cur);
    }

    return std::make_tuple(name, seq);
}


int main() {
    // verify input
    std::string filename("/mnt/d/Homo_sapiens.GRCh38.dna.chromosome.1.fa/Homo_sapiens.GRCh38.dna.chromosome.1.fa");
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return 1;
    }

    // std::vector<std::map<char, int>> counts; // counts[i][c] is the number of times character c appears in column i
    std::string reference_seq;

    // read allignments and print them
    std::string name, seq;
    while (!in.eof()) {
        std::tie(name, seq) = next_sample(in);

        // update columns and counts size
        if (columns == 0) {
            columns = seq.size();
        }

    }

    // determain column i is parsimony informative
    std::vector<bool> informative(counts.size(), false);

    // find number of characters that occurs at least twice
    for (int i = 0; i < columns; ++i) {
        int n = 0;
        for (int c = 0; c < counts.size(); ++c) {
            if (counts[c][i] > 1) {
                ++n;
            }
        }
        if (n > 1) {
            informative[i] = true;
        }
    }

    // print informative columns
    for (int i = 0; i < counts.size(); ++i) {
        if (informative[i]) {
            std::cout << i << std::endl;
        }
    }
}