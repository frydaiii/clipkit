#include <tuple>
#include <fstream>
#include "const.h"

int is_unknown_site(char base, std::string seq_type);
std::tuple<std::string, std::string> next_sample(std::ifstream &in);