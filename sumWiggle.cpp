#include "wiggle_tools.hpp"
#include "bed_tools.hpp"
#include "tools.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
    std::map <char, std::string> args = getArgs(std::vector<std::string>(argv, argv + argc));
    std::string wiggle_filename, AOE_filename, output_filename;
    try {
        wiggle_filename = args.at('w');
        AOE_filename = args.at('a');
        output_filename = args.at('o');
    } catch (std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are : \n";
        std::cout << "\t+a AOE file \n";
        std::cout << "\t+w wiggle file \n";
        std::cout << "\t+o output file \n";
        throw;
    }

    std::cout << "Reading files" << std::endl;
    AOE_file intervals(AOE_filename, read);
    intervals.readWholeFile();
    wiggle_file values(wiggle_filename, read);
    values.readWholeFile();
    std::ofstream output_file(output_filename);

    std::map <int, std::vector<double>> summed_values;

    std::cout << "Intersecting and counting" << std::endl;
    for(const auto& result: intervals.intersect(values)) {
        std::vector <valid_double> values(dynamic_cast <const wiggle_entry*>(result.hit) -> getSubset(&result.result));
        for(int i(0); i < values.size(); i++) {
            if(summed_values[dynamic_cast <AOE_entry*> (result.source) -> getRelativePos(result.result.getStart() + i)].size() < 2) {
                summed_values[dynamic_cast <AOE_entry*> (result.source) -> getRelativePos(result.result.getStart() + i)] = std::vector <double>(2, 0);
            }
            if(values[i].is_valid) {
                if(values[i].value > 1.0) {
                    std::cout << "before : " << values[i - 1].value << std::endl;
                    std::cout << "error is a line " << i << std::endl;
                    std::cout << "value is : " << values[i].value << std::endl;
                    std::cout << result.hit->getString() << std::endl;
                    throw(666);
                }
                summed_values[dynamic_cast <AOE_entry*> (result.source) -> getRelativePos(result.result.getStart() + i)][0] += values[i].value;
                summed_values[dynamic_cast <AOE_entry*> (result.source) -> getRelativePos(result.result.getStart() + i)][1] += 1.0;
            }
        }
    }

    std::cout << "Writing results" << std::endl;
    for(const auto& pair: summed_values) {
        output_file << pair.first << '\t' << pair.second[0] << '\t' << int(pair.second[1]) << '\n';
    }
    return 0;
}