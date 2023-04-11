#include "snp-sites.h"

std::tuple<std::string, std::vector<int>, std::vector<int>> snp_sites(std::string input_file) {
    std::string seq_type = "";
    // verify input
    std::ifstream in(input_file);
    if (!in) {
        std::cerr << "Error: cannot open file " << input_file << std::endl;
        return std::make_tuple(seq_type, std::vector<int>(), std::vector<int>());
    }


    // init reference_seq and gap_nums
    int num_of_seqs = 1; // start at 1 because we already read the first sequence
    std::string name, seq;
    std::tie(name, seq) = next_sample(in);
    std::string reference_seq = seq;
    std::vector<int> gap_nums(seq.size(), 0);

    for (int i = 0; i < seq.size(); i++) {
        if (is_unknown_site(seq[i], seq_type)) {
            gap_nums[i]++;
        }
    }


    std::unordered_set<char> sites; // sites is a set of all bases that occur at least once
    while (!in.eof()) {
        std::tie(name, seq) = next_sample(in);
        num_of_seqs++;

        // if reference_seq haven't been initialized
        if (reference_seq.size() == 0) {
            // init reference_seq with seq's length of N characters
            reference_seq = std::string(seq.size(), 'N');
        }

        // find snp-sites
        for (int i = 0; i < seq.size(); i++) {
            // ignore snp-sites and unknown sites
            if (reference_seq[i] == '>') continue;
            if (is_unknown_site(seq[i], seq_type)) {
                gap_nums[i]++;
                continue;
            }

            // update sites
            sites.insert(seq[i]);

            // check if site i is a snp-site
            if (!is_unknown_site(reference_seq[i], seq_type) && reference_seq[i] != seq[i]) {
                reference_seq[i] = '>';
            } else if (is_unknown_site(reference_seq[i], seq_type) && reference_seq[i] != seq[i]) {
                reference_seq[i] = seq[i];
            }
        }
    }

    // get snp-site and constant-sites locations
    std::vector<int> snp_sites_loc;
    std::vector<int> constant_sites_loc;
    for (int i = 0; i < reference_seq.size(); i++) {
        if (reference_seq[i] == '>') {
            snp_sites_loc.push_back(i);
        } else {
            if (gap_nums[i] < num_of_seqs - 1) { // at least one base occurs at least twice
                constant_sites_loc.push_back(i);
            }
        }
    }

    in.close();

    // determine sequence type
    if (sites.size() > 5) {
        seq_type = AA;
    } else {
        seq_type = NT;
    }

    return std::make_tuple(seq_type, snp_sites_loc, constant_sites_loc);
}