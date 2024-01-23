#include "bed_tools.hpp"
#include "tools.hpp"
#include <iostream>

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
        std::cout << "\t+s count values by strand (hit & result)\n";
        std::cout << "\t+p start/stop/mid/total\n";
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

    bool count_ids(false);
    try {
        args.at('d');
        count_ids = true;
    } catch(std::out_of_range) {
        std::cout << "counting for whole bed" << std::endl;
    }

    bool count_strand(false);
    try {
        args.at('s');
        count_strand = true;
    } catch(std::out_of_range) {
        std::cout << "won't keep strand" << std::endl;
    }

    enum type_count {whole, start, stop, mid};
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
        }
    } catch(std::out_of_range) {
        std::cout << "Counting for whole interval" << std::endl;
    }
    
    std::cout << "Loading AOE" << std::endl;
    AOE_file AOEs(AOE_filename, read);
    AOEs.readWholeFile();

    std::ofstream output_file(output_filename);
    std::map <std::string, std::map <char, std::map <char, std::map<int, int>>>> summed_values;

    std::cout << "Loading bed" << std::endl;
    bed_file ints_to_count(bed_filename, read);
    
    std::cout << "Intersecting" << std::endl;
    std::vector <intersect_results> results;
    while(ints_to_count.remainToRead()) {
        if(memory_saving) {
            ints_to_count.eraseAndLoadBlock(number_blocks);
        } else {
            ints_to_count.readWholeFile();
        }
        results = AOEs.intersect(ints_to_count);
        for(const auto& entry: results) {
            std::string id("no_id");
            char strand_source('+');
            char strand_hit('+');
            if(count_ids) {
                id = entry.result.getID();
            }
            if(count_strand) {
                strand_source = entry.source -> getStrand();
                strand_hit = entry.hit -> getStrand();
            }
            if(count == type_count::whole) {
                for(int i(entry.result.getStart()); i < entry.result.getEnd(); i++) {
                    summed_values[id][strand_source][strand_hit][dynamic_cast<AOE_entry*>(entry.source) -> getRelativePos(i)] ++;
                }
            } else if(count == type_count::start) {
                summed_values[id][strand_source][strand_hit][dynamic_cast<AOE_entry*>(entry.source) -> getRelativePos(entry.hit->getStart())] ++;
            } else if(count == type_count::stop) {
                summed_values[id][strand_source][strand_hit][dynamic_cast<AOE_entry*>(entry.source) -> getRelativePos(entry.hit->getEnd())] ++;
            } else if(count == type_count::mid) {
                summed_values[id][strand_source][strand_hit][dynamic_cast<AOE_entry*>(entry.source) -> getRelativePos((entry.hit->getStart() + entry.hit->getEnd())/2)] ++;
            }
        }
    }    

    std::cout << "Writing results" << std::endl;
    std::string tmp("");
    for(const auto &idToMap: summed_values) { // pair id:map
        for(const auto &strandSourceToMap: idToMap.second) { // pair strand source : map
            for(const auto &strandHitToMap: strandSourceToMap.second) {
                for(const auto &posToInt: strandHitToMap.second) {
                    tmp = "";
                    if(count_ids) {
                        tmp += idToMap.first + "\t";
                    }
                    if(count_strand) {
                        tmp += std::string(1, strandSourceToMap.first) + "\t" + std::string(1, strandHitToMap.first) + "\t";
                    }
                    tmp += std::to_string(posToInt.first) + "\t" + std::to_string(posToInt.second) + "\n";
                    output_file << tmp;
                }
            }
        }
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