#include <tuple>
#include <fstream>
#include "const.h"

int is_unknown_site(char, std::string);
std::tuple<std::string, std::string> next_sample(std::ifstream &);