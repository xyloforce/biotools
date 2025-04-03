#include "bed_tools.hpp"
#include "tools.hpp"
#include <iostream>

// FIXME implement correctly saving by chromosome : change FILENAME

void write_results(std::string output_filename, std::map <std::string, std::map <int, int>> &values_map) {
    std::ofstream output_file(output_filename);
    for (const auto& idToMap : values_map) { // pair id:map
        for (const auto& posToMap : idToMap.second) { // pair strand source : map
            if (idToMap.first != "") {
				output_file << idToMap.first;
            }
            output_file << std::to_string(posToMap.first) + "\t" + std::to_string(posToMap.second) + "\n";
			// id : pos : value
        }
    }
    values_map.clear(); // empty the container bc contents have been written
}

int main(int argc, char *argv[]) {
    std::map <char, std::string> args = getArgs(std::vector<std::string>(argv, argv + argc));
    std::string bed_filename, AOE_filename, output_filename;
    try {
        bed_filename = args.at('b');
        AOE_filename = args.at('a');
        output_filename = args.at('o');
    } catch (std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are : \n";
        std::cout << "\t+a AOE file \n";
        std::cout << "\t+b bed file \n";
        std::cout << "\t+o output file \n";
        std::cout << "Optionnal parameters are :\n";
        std::cout << "\t+m save memory but slower\n";
        std::cout << "\t+i number of lines to load at a time\n";
        std::cout << "\t+d count values by bed identifier\n";
        std::cout << "\t+k id : keep both, source or hit" << std::endl;
        std::cout << "\t+s count values by strand (hit & result)\n";
		std::cout << "\t+v only match same strand\n";
        std::cout << "\t+p start/stop/mid/whole\n";
        std::cout << "\t+c save results by chromosome\n";
        throw;
    }
    
    bool memory_saving(false);
    int number_blocks(100000);
    try {
        args.at('m');
        memory_saving = true;
        try {
            number_blocks = stoi(args.at('i'));
        } catch(std::out_of_range) {
            std::cout << "using default number of blocks" << std::endl;
        }
    } catch(std::out_of_range) {
        std::cout << "memory saving not set" << std::endl;
    }

    id_status status = both;
    bool count_ids(false);
    try {
        args.at('d');
        count_ids = true;
        try {
            std::string statusS(args.at('k'));
            if (statusS == "source") {
                status = source;
            }
            else if (statusS == "hit") {
                status = hit;
            }
            else {
                status = both;
            }
        }
        catch (std::out_of_range) {
            status = both;
        }
    }
    catch (std::out_of_range) {
        std::cout << "Ignoring ids" << std::endl;
    }

    bool count_by_chromosome(false);
    try {
        args.at('c');
        count_by_chromosome = true;
    }
    catch (std::out_of_range) {
        // do nothing bc it means that arg is not set
    }

    bool count_strand(false);
    try {
        args.at('s');
        count_strand = true;
    } catch(std::out_of_range) {
        std::cout << "won't keep strand" << std::endl;
    }
	
	bool only_strand(false);
    try {
        args.at('v');
        only_strand = true;
    } catch(std::out_of_range) {
        std::cout << "won't keep strand" << std::endl;
    }

    enum type_count {whole, start, stop, mid}; // what shall we count : full interval ? only start / stop ? Only the middle ?
    type_count count(type_count::whole);
    try {
        std::string val(args.at('p'));
        if(val == "whole") {
            count = type_count::whole;
        } else if(val == "start") {
            count = type_count::start;
        } else if(val == "stop") {
            count = type_count::stop;
        } else if(val == "mid") {
            count = type_count::mid;
        } else {
            std::cout << val << " : invalid option" << std::endl;
            throw std::invalid_argument("change p parameter");
        }
    } catch(std::out_of_range) {
        std::cout << "Counting for whole interval" << std::endl;
    }
    
    std::cout << "Loading AOE" << std::endl;
    AOE_file AOEs(AOE_filename, read);
    AOEs.readWholeFile();

    //std::ofstream output_file(output_filename);
    //std::map <std::string, std::map <char, std::map <char, std::map<int, int>>>> summed_values;
    std::map <std::string, std::map <int, int>> summed_values;

    std::cout << "Loading bed" << std::endl;
    bed_file ints_to_count(bed_filename, read);
    
    std::cout << "Intersecting" << std::endl;
    std::vector <intersect_results> results;
    std::string chromosome("");
    while(ints_to_count.remainToRead()) {
        if(memory_saving) {
            ints_to_count.eraseAndLoadBlock(number_blocks);
        } else {
            ints_to_count.readWholeFile();
        }
        results = AOEs.intersect(ints_to_count, only_strand, status);
        
        for(const auto& entry: results) {
            if (chromosome != "" && count_by_chromosome && chromosome != entry.source -> getChr()) {
                std::string filename(output_filename + "_" + chromosome + ".tsv");
                write_results(filename, summed_values);
            }
            chromosome = entry.source -> getChr();
            std::string key("");
            /*std::string id("no_id");
            char strand_source('+');
            char strand_hit('+');*/
            if(count_ids) {
                key += entry.result.getID();
                if(key.find("._") == 0) {
                    key = key.erase(0, 2);
                }
				if(key == ".") { // if i asked for an id it would be nice to have a real one
                    key = ""; // delete the point bc we don't need it
                    if(status == source) {
                        key += entry.source -> getChr() + std::to_string(entry.source -> getStart()) + std::to_string(entry.source -> getEnd());
                    } else if(status == hit) {
                        key += entry.hit -> getChr() + std::to_string(entry.hit -> getStart()) + std::to_string(entry.hit -> getEnd());
                    } else {
                        key += entry.result.getChr() + std::to_string(entry.result.getStart()) + std::to_string(entry.result.getEnd());
                    }
                }
				key += "\t";
            }
            if(count_strand) {
                key += entry.source -> getStrand();
                key += "\t";
                key += entry.hit -> getStrand();
				key += "\t";
            }

            if(count == type_count::whole) {
                for(int i(entry.result.getStart()); i < entry.result.getEnd(); i++) {
                    summed_values[key][dynamic_cast<AOE_entry*>(entry.source) -> getRelativePos(i)] ++;
                }
            } else if(count == type_count::start) {
                if(entry.hit -> getStrand() == '+') {
                    summed_values[key][dynamic_cast<AOE_entry*>(entry.source) -> getRelativePos(entry.hit->getStart())] ++;
                } else {
                    summed_values[key][dynamic_cast<AOE_entry*>(entry.source) -> getRelativePos(entry.hit->getEnd())] ++;
                }
            } else if(count == type_count::stop) {
                if(entry.hit -> getStrand() == '-') {
                    summed_values[key][dynamic_cast<AOE_entry*>(entry.source) -> getRelativePos(entry.hit->getStart())] ++;
                } else {
                    summed_values[key][dynamic_cast<AOE_entry*>(entry.source) -> getRelativePos(entry.hit->getEnd())] ++;
                }
            } else if(count == type_count::mid) {
                summed_values[key][dynamic_cast<AOE_entry*>(entry.source) -> getRelativePos((entry.hit->getStart() + entry.hit->getEnd())/2)] ++;
            }
        }
    }    
    if (count_by_chromosome) {
        std::string filename(output_filename + "_" + chromosome + ".tsv");
        std::cout << filename << std::endl;
        write_results(filename, summed_values);
    }
    else {
        std::cout << "Writing results" << std::endl;
        write_results(output_filename, summed_values);
    }
    return 0;
}

/*
objective : get back AOEs from the intersection bc they have a zero
then its only i from start - zero to end - zero
however intersection can be bio_entries
but should return either bio_entry or AOE_entries
or :
- create array of bio_entry unique ptr pointing towards AOE_entries
- intersect
- return AOE_entries
*/
