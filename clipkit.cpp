#include "clipkit.h"

std::vector<std::map<char, int>> num_ocurrences_snp(std::string input_file, std::string seq_type, std::vector<int> snp_sites_loc) {
    // open snp-sites file for reading
    std::ifstream in(input_file);

    std::vector<std::map<char, int>> base_counts(snp_sites_loc.size()); // base_counts[i] is a map of base counts at snp_sites_loc[i] of the sequence.
    std::string name, seq;
    while (!in.eof()) {
        std::tie(name, seq) = next_sample(in);

        for (int i = 0; i < snp_sites_loc.size(); i++) { // iterate through snp-sites
            int j = snp_sites_loc[i];
            // process site 
            if (is_unknown_site(seq[j], seq_type)) {
                // ignore unknown sites
                continue;
            } else {
                // update base_counts[i]
                if (base_counts[i].find(seq[j]) == base_counts[i].end()) {
                    // new base
                    base_counts[i][seq[j]] = 1;
                } else {
                    // existing base
                    base_counts[i][seq[j]]++;
                }
            }
        }
    }

    return base_counts;
}

std::vector<int> determine_parsimony_informative(std::vector<int> snp_sites_loc, std::vector<std::map<char, int>> base_counts) {
    std::vector<int> parsimony_informative_sites_loc;

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
            parsimony_informative_sites_loc.push_back(snp_sites_loc[i]);
        }
    }

    return parsimony_informative_sites_loc;
}

void write_result(std::string input_file, std::vector<int> sites_loc, std::string output) {
    // open snp-sites file for reading
    std::ifstream in(input_file);
    
    // print parsimony-informative sites to file
    std::ofstream out(output);
    std::string name, seq;
    while (!in.eof()) {
        std::tie(name, seq) = next_sample(in);

        // print sample name
        out << ">" << name << std::endl;

        // print parsimony-informative sites, newline at 60th site
        int count = 0;
        for (int i = 0; i < sites_loc.size(); i++) {
            count++;
            if (count % 60 == 1 && count > 1) {
                out << std::endl;
            }
            out << seq[sites_loc[i]];
        }
        if (count > 0) out << std::endl;
    }

    // close files
    in.close();
    out.close();
}

void clipkit(std::string input_file, std::string output_file) {
    // get sequence type, snp-sites locations and constant-sites locations
    std::string seq_type;
    std::vector<int> snp_sites_loc;
    std::vector<int> const_sites_loc;
    std::tie(seq_type, snp_sites_loc, const_sites_loc) = snp_sites(input_file);

    // find parsimony-informative sites
    std::vector<std::map<char, int>> base_snp_counts = num_ocurrences_snp(input_file, seq_type, snp_sites_loc); 
    std::vector<int> pi_sites_loc = determine_parsimony_informative(snp_sites_loc, base_snp_counts);

    // print result
    // combine parsimony-informative sites and constant sites into new vector
    std::vector<int> sites_loc(pi_sites_loc.size() + const_sites_loc.size());
    std::merge(pi_sites_loc.begin(), pi_sites_loc.end(), const_sites_loc.begin(), const_sites_loc.end(), sites_loc.begin());
    std::sort(sites_loc.begin(), sites_loc.end());

    write_result(input_file, pi_sites_loc, output_file);
}
