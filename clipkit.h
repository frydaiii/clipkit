#include "snp-sites.h"
#include <map>

/*
    Return number of occurences of each base.
*/
std::vector<std::map<char, int>> num_ocurrences(std::string input, std::string seq_type);

/*
    Return true if site is parsimony-informative.
*/
std::vector<bool> determine_parsimony_informative(std::vector<std::map<char, int>> base_counts);

/*
    Write result to file.
*/
void write_result(std::string snp_file, std::vector<bool> is_parsimony_informative_site, std::string output);

void clipkit(std::string input_filename, std::string output_filename);