#include "snp-sites.h"
#include <map>

/*
    Return number of occurences of each snp-site base.
*/
std::vector<std::map<char, int>> num_ocurrences_snp(std::string, std::string, std::vector<int>);

/*
    Return locations parsimony-informative sites.
*/
std::vector<int> determine_parsimony_informative(std::vector<int>, std::vector<std::map<char, int>>);

/*
    Write result to file.
*/
void write_result(std::string, std::vector<int>, std::string);

void clipkit(std::string, std::string);