#include "vcf_tools.hpp"
#include "bed_tools.hpp"
#include "fasta_tools.cpp"
#include "tools.hpp"

int main(int argc, char* argv[]) {
    std::map <char, std::string> args = getArgs(std::vector <std::string> (argv, argv + argc));
    std::string vcf_filename, fasta_input_filename, vcf_output_filename;
    std::cout << "Starting" << std::endl;
    try {
        vcf_filename = args.at('v');
        fasta_input_filename = args.at('f');
        vcf_output_filename = args.at('o');
    } catch (std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are :\n";
        std::cout << "\t+v vcf filename\n";
        std::cout << "\t+f input fasta\n";
        std::cout << "\t+o output vcf\n";
        throw std::out_of_range("Missing arguments");
    }

    std::cout << "Loading files" << std::endl;

    vcf_file mutations(vcf_filename, read), output_file(vcf_output_filename, write);
    fasta_file input_fasta(fasta_input_filename, read, standard);

    input_fasta.readWholeFile();

    std::cout << "Checking mutations" << std::endl;
    bool ok(true);
    vcf_entry* vcf_converted(0);
    std::unique_ptr <bio_entry> entry;
    int valid(0), invalid(0), out_of_bounds(0), non_attributed(0);

    while(mutations.remainToRead()) {
        ok = true;
        try {
            entry = std::move(mutations.readLine());
            vcf_converted = dynamic_cast <vcf_entry*>(entry.get());
            if(!vcf_converted) {
                throw std::bad_cast();
            }
        } catch(std::out_of_range) {
            ok = false;
            std::cout << "skipped a empty entry, probably whiteline at the eof" << std::endl;
        }
        if(ok) {
            try {
                fasta_entry* fasta_converted = dynamic_cast <fasta_entry*> (input_fasta.getEntriesByChr(vcf_converted -> getChr())[0]);
                try {
                    if(fasta_converted -> getBase(vcf_converted -> getStart()) == vcf_converted -> getRef()[0]) {
                        valid ++;
                        output_file.writeBioLine(*vcf_converted);
                    } else {
                        invalid ++;
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