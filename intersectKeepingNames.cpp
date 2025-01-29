#include "bed_tools.hpp"
#include "vcf_tools.hpp"
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
        std::cout << "Optionnal :\n";
        std::cout << "\t+m id : keep both, source or hit\n";
        std::cout << "\t+v masked vcf\n";
        throw;
    }

    id_status status = source;
    try {
        std::string statusS(args.at('m'));
        if(statusS == "source") {
            status = source;
        } else if(statusS == "hit") {
            status = hit;
        } else {
            status = both;
        }
    } catch (std::out_of_range) {
        status = source;
    }

    bool is_a_vcf(false);
    try  {
        args.at('v');
        is_a_vcf = true;
    } catch (std::out_of_range) {
        // nothing bc then default is bed
    }

    std::cout << "Loading beds..." << std::endl;
    bed_file id_source(bed_filename, read);
    id_source.readWholeFile();

    bio_file *mask;
    vcf_file vcf_temp(mask_filename, none);
    bed_file bed_temp(mask_filename, none);
    if (is_a_vcf) {
        mask = &vcf_temp;
    } else {
        mask = &bed_temp;
    }
    mask -> typeToRead("");
    mask -> readWholeFile();
    std::cout << "Intersecting" << std::endl;
    id_source.apply_intersect(*mask, false, status);
    std::cout << "Write results" << std::endl;
    id_source.typeToWrite(output_filename);
    id_source.writeToFile();
}
