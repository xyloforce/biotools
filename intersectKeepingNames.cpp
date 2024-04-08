#include "bed_tools.hpp"
#include "tools.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
    std::map <char, std::string> args = getArgs(std::vector<std::string>(argv, argv + argc));
    std::string bed_filename, mask_filename, output_filename;
    try {
        bed_filename = args.at('a');
        mask_filename = args.at('b');
        output_filename = args.at('o');
    } catch (std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are : \n";
        std::cout << "\t+a bed file \n";
        std::cout << "\t+b mask bed file \n";
        std::cout << "\t+o output file \n";
        throw;
    }
    std::cout << "Loading beds..." << std::endl;
    bed_file id_source(bed_filename, read);
    id_source.readWholeFile();
    bed_file mask(mask_filename, read);
    mask.readWholeFile();
    std::cout << "Intersecting" << std::endl;
    id_source.apply_intersect(mask);
    std::cout << "Write results" << std::endl;
    id_source.typeToWrite(output_filename);
    id_source.writeToFile();
}