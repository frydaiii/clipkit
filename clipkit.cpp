#include <string>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>
#include <map>
#include <unordered_set>
#include <filesystem>
namespace fs = std::filesystem;

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

const std::string AA = "aa", NT = "nt";

int is_unknown_site(char base, std::string seq_type)
{
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

void clipkit(std::string input_filename, std::string output_filename) {
    std::string seq_type = "";
    // verify input
    std::ifstream in(input_filename);
    if (!in) {
        std::cerr << "Error: cannot open file " << input_filename << std::endl;
        return;
    }

    // read allignments and find snp-sites
    std::string name, seq;
    std::string reference_seq;
    std::unordered_set<char> sites;

    while (!in.eof()) {
        std::tie(name, seq) = next_sample(in);

        // if reference_seq haven't been initialized
        if (reference_seq.size() == 0) {
            // init reference_seq with seq's length of N characters
            reference_seq = std::string(seq.size(), 'N');
        }

        // find snp-sites
        for (int i = 0; i < seq.size(); i++) {
            // ignore snp-sites and unknown sites
            if (reference_seq[i] == '>' || is_unknown_site(seq[i], seq_type)) continue;

            // update sites
            sites.insert(seq[i]);

            // check if site i is a snp-site
            if (reference_seq[i] != seq[i]) {
                reference_seq[i] = '>';
            }
        }
    }

    // determine sequence type
    if (sites.size() > 5) {
        seq_type = AA;
    } else {
        seq_type = NT;
    }

    // get snp-site locations
    std::vector<int> snp_sites_loc;
    for (int i = 0; i < reference_seq.size(); i++) {
        if (reference_seq[i] == '>') {
            snp_sites_loc.push_back(i);
        }
    }

    // reset file pointer
    in.clear();
    in.seekg(0, std::ios_base::beg);

    // print snp-sites to file
    std::ofstream out("snp_sites.fa");
    while (!in.eof()) {
        std::tie(name, seq) = next_sample(in);

        // print sample name
        out << ">" << name << std::endl;
        
        // print snp-sites
        for (int i = 0; i < snp_sites_loc.size(); i++) {
            out << seq[snp_sites_loc[i]];
        }
        out << std::endl;
    }

    // close original input file
    in.close();
    out.close();
    
    // open snp-sites file for reading
    in.open("snp_sites.fa");

    // find parsimony-informative sites
    std::vector<bool> is_parsimony_informative_site(snp_sites_loc.size(), false); // is_parsimony_informative_site[i] = true if seq[i] is a parsimony-informative site
    std::vector<std::map<char, int>> base_counts(snp_sites_loc.size()); // base_counts[i] is a map of base counts at snp_sites_loc[i]
    while (!in.eof()) {
        std::tie(name, seq) = next_sample(in);

        // find parsimony-informative sites
        for (int i = 0; i < seq.size(); i++) {
            // process site i
            if (is_unknown_site(seq[i], seq_type)) {
                // ignore unknown sites
                continue;
            } else {
                // update base_counts[i]
                if (base_counts[i].find(seq[i]) == base_counts[i].end()) {
                    // new base
                    base_counts[i][seq[i]] = 1;
                } else {
                    // existing base
                    base_counts[i][seq[i]]++;
                }
            }
        }
    }
    for (int i = 0; i < base_counts.size(); i++) {
        // counting number of bases that occur at least twice
        int count = 0;
        for (auto it = base_counts[i].begin(); it != base_counts[i].end(); it++) {
            if (it->second >= 2) {
                count++;
            }
        }
        // if there are at least two bases that occur at least twice
        if (count >= 2) {
            is_parsimony_informative_site[i] = true;
        }
    }

    //reset file pointer
    in.clear();
    in.seekg(0, std::ios_base::beg);

    // print parsimony-informative sites to file
    out.open(output_filename);
    while (!in.eof()) {
        std::tie(name, seq) = next_sample(in);

        // print sample name
        out << ">" << name << std::endl;

        // print parsimony-informative sites, newline at 60th site
        int count = 0;
        for (int i = 0; i < is_parsimony_informative_site.size(); i++) {
            if (is_parsimony_informative_site[i] == true) {
                count++;
                if (count % 60 == 1 && count > 1) {
                    out << std::endl;
                }
                out << seq[i];
            }
        }
        if (count > 0) out << std::endl;
    }

    // close files
    in.close();
    out.close();
}

int main() {
    // get all files in current directory
    std::string path = "/mnt/d/GBE_Shen_etal_2016/GBE_2016/Full_length_dataset/Alignments_mammals/aa";
    for (const auto & entry_it : fs::directory_iterator(path))
        if (entry_it.is_regular_file()) {
            std::string input_filename = entry_it.path().string();

            // gen output file by change input file's extension
            std::string extension = fs::path(input_filename).extension().string();

            // ignore non-fasta files
            if (extension != ".fasta") continue;
            std::string new_extension = ".fasta.clipkit";
            std::string output_filename = input_filename.substr(0, input_filename.size() - extension.size()) + new_extension;

            // run clipkit
            clipkit(input_filename, output_filename);
        }
    return 0;
}