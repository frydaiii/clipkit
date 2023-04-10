#include "clipkit.h"

std::vector<std::map<char, int>> num_ocurrences(std::string input, std::string seq_type) {
    // open snp-sites file for reading
    std::ifstream in(input);

    // find parsimony-informative sites
    std::vector<std::map<char, int>> base_counts; // base_counts[i] is a map of base counts at snp_sites_loc[i]
    std::string name, seq;
    while (!in.eof()) {
        std::tie(name, seq) = next_sample(in);

        // init is_parsimony_informative_site and base_counts
        if (base_counts.size() == 0) {
            base_counts = std::vector<std::map<char, int>>(seq.size());
        }

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

    return base_counts;
}

std::vector<bool> determine_parsimony_informative(std::vector<std::map<char, int>> base_counts) {
    std::vector<bool> is_parsimony_informative_site(base_counts.size(), false); // is_parsimony_informative_site[i] = true if seq[i] is a parsimony-informative site

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

    return is_parsimony_informative_site;
}

void write_result(std::string snp_file, std::vector<bool> is_parsimony_informative_site, std::string output) {
    // open snp-sites file for reading
    std::ifstream in(snp_file);
    
    // print parsimony-informative sites to file
    std::ofstream out(output);
    std::string name, seq;
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

void clipkit(std::string input_file, std::string output_file) {
    std::string snp_file = "snp_sites.fa"; // temporary file for snp-sites output
    std::string seq_type = snp_sites(input_file, snp_file);

    // find parsimony-informative sites
    std::vector<std::map<char, int>> base_counts = num_ocurrences(snp_file, seq_type); 
    std::vector<bool> pi_sites = determine_parsimony_informative(base_counts);

    // print result
    write_result(snp_file, pi_sites, output_file);

    // remove snp_sites.fa
    remove(snp_file.c_str());
}
