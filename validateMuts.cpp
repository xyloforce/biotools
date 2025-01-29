#include "vcf_tools.hpp"
#include "bed_tools.hpp"
#include "fasta_tools.cpp"
#include "tools.hpp"

int main(int argc, char* argv[]) {
    std::map <char, std::string> args = getArgs(std::vector <std::string> (argv, argv + argc));
    std::string vcf_filename, fasta_input_filename;
    std::cout << "Starting" << std::endl;
    try {
        vcf_filename = args.at('v');
        fasta_input_filename = args.at('f');
    } catch (std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are :\n";
        std::cout << "\t+v vcf filename\n";
        std::cout << "\t+f input fasta\n";
        std::cout << "Optionnal :" << std::endl;
        std::cout << "\t+i consider the alt instead of the ref base" << std::endl;
        std::cout << "\t+1 validated mutations vcf\n";
        std::cout << "\t+2 non-validated mutations vcf\n";
        throw std::out_of_range("Missing arguments");
    }

    std::cout << "Loading files" << std::endl;

    std::string validated_filename, invalidated_filename;
    vcf_file mutations(vcf_filename, read), output_validated, output_invalidated;
    fasta_file input_fasta(fasta_input_filename, read, standard);

    input_fasta.readWholeFile();

    std::cout << "Checking mutations" << std::endl;
    bool ok(true);
    vcf_entry* vcf_converted(0);
    std::unique_ptr <bio_entry> entry;
    int valid(0), invalid(0), out_of_bounds(0), non_attributed(0);

    bool consider_alt(false);
    try {
        args.at('i');
        std::cout << "Considering ALT base" << std::endl;
        consider_alt = true;
    } catch(std::out_of_range) {
        std::cout << "Considering REF base" << std::endl;
    }

    enum type_keep {validated, invalidated, none, both};
    type_keep keep(none);
    try {
        validated_filename = args.at('1');
        keep = validated;
        output_validated.typeToWrite(validated_filename);
    } catch(std::out_of_range) {
        std::cout << "Not keeping validated mutations" << std::endl;
    }
    try {
        invalidated_filename = args.at('2');
        output_invalidated.typeToWrite(invalidated_filename);
        if(keep == validated) {
            std::cout << "Keeping both kinds" << std::endl;
            keep = both;
        } else {
            keep = invalidated;
        }
    } catch(std::out_of_range) {
        std::cout << "Not keeping invalidated mutations" << std::endl;
    }

    while(mutations.remainToRead()) {
        ok = true;
        try {
            entry = std::move(mutations.readLine());
            vcf_converted = dynamic_cast <vcf_entry*>(entry.get());
            if(!vcf_converted) {
                throw std::bad_cast();
            }
        } catch(const std::out_of_range& e) {
            ok = false;
            std::cout << "invalid line: " << e.what() << std::endl;
        }
        if(ok) {
            try {
                fasta_entry* fasta_converted = dynamic_cast <fasta_entry*> (input_fasta.getEntriesByChr(vcf_converted -> getChr())[0]);
                char current_base = vcf_converted -> getRef()[0];
                if(consider_alt) {
                    current_base = vcf_converted -> getAlt(0)[0];
                }
                try {
                    if(fasta_converted -> getBase(vcf_converted -> getStart()) == current_base) {
                        valid ++;
                        if(keep == validated || keep == both) {
                            output_validated.writeBioLine(*vcf_converted);
                        }
                    } else {
                        invalid ++;
                        if(keep == invalidated || keep == both) {
                            output_invalidated.writeBioLine(*vcf_converted);
                        }
                    }
                } catch(std::out_of_range) {
                    out_of_bounds ++;
                }
            } catch(std::out_of_range) {
                non_attributed ++;
            }
        }
    }
    std::cout << "Valid : " << valid << std::endl;
    std::cout << "Invalid : " << invalid << std::endl;
    std::cout << "Out of bounds : " << out_of_bounds << std::endl;
    std::cout << "Chromosome non present in the genome : " << non_attributed << std::endl;
}