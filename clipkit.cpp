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

int is_unknown_site(char base)
{
    switch (base) {
        case 'N':
        case 'n':
        case '-':
        case '?':
            return 1;
        default:
            return 0;
    }
}

int main() {
    // verify input
    std::string filename("/mnt/d/Homo_sapiens.GRCh38.dna.chromosome.1.fa/Homo_sapiens.GRCh38.dna.chromosome.1.fa");
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return 1;
    }


    // read allignments and find snp-sites
    std::string name, seq;
    std::string reference_seq;

    while (!in.eof()) {
        std::tie(name, seq) = next_sample(in);

        // if reference_seq haven't been initialized
        if (reference_seq.size() == 0) {
            // init reference_seq with seq's length of N characters
            reference_seq = std::string(seq.size(), 'N');
        }

        // find snp-sites
        for (int i = 0; i < reference_seq.size(); i++) {
            // ignore snp-sites
            if (reference_seq[i] == '>') continue;

            // process site i
            if (reference_seq[i] == 'N') {
                // update reference_seq[i] if needed
                if (!is_unknown_site(seq[i])) {
                    reference_seq[i] = seq[i];
                }
            } else {
                // check if site i is a snp-site
                if (!is_unknown_site(seq[i]) && reference_seq[i] != seq[i]) {
                    reference_seq[i] = '>';
                }
            }
        }
    }

    // reset file pointer
    in.clear();
    in.seekg(0, std::ios_base::beg);

    // write snp-sites to file
    std::ofstream out("snp-sites.fa");
    while (!in.eof()) {
        std::tie(name, seq) = next_sample(in);

        out << ">" << name << std::endl;
        for (int i = 0; i < reference_seq.size(); i++) {
            if (reference_seq[i] == '>') {
                out << seq[i];
            }
        }
        out << std::endl;
    }

    // close files
    in.close();
    out.close();
}