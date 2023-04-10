#include "clipkit.h"
#include <filesystem>

namespace fs = std::filesystem;

int main() {
    // get all files in current directory
    // std::string path = "/mnt/d/GBE_Shen_etal_2016/GBE_2016/Full_length_dataset/Alignments_mammals/aa";
    std::string path = "/mnt/d/Homo_sapiens.GRCh38.dna.chromosome.1.fa";
    for (const auto & entry_it : fs::directory_iterator(path))
        if (entry_it.is_regular_file()) {
            std::string input_file = entry_it.path().string();

            // gen output file by change input file's extension
            std::string extension = fs::path(input_file).extension().string();

            // ignore non-fasta files
            if (extension != ".fasta" && extension != ".fa") continue;
            std::string new_extension = ".fasta.clipkit";
            std::string output_file = input_file.substr(0, input_file.size() - extension.size()) + new_extension;

            // run clipkit
            clipkit(input_file, output_file);
        }
    return 0;
}