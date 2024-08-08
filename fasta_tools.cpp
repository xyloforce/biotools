#include "fasta_tools.hpp"
#include <iostream>
#include <string_view>
#include <sstream>
#include <regex>

fasta_entry::fasta_entry(std::string chr, long start, long end, std::string id, char strand, std::string sequence): bio_entry(chr, start, end, strand, id), m_sequence(sequence) {
    if(end - start != m_sequence.size()) {
        std::cout << "sequence is : " << sequence.size() << " bp\n";
        std::cout << "end - start is " << end - start << " bp\n";
        throw std::invalid_argument("sequence is not the same size as its coordinates");
    }
}

fasta_entry::fasta_entry(bio_entry entry, std::string sequence): bio_entry(entry), m_sequence(sequence) {

}

fasta_entry::fasta_entry(): bio_entry(), m_sequence("") {

}

fasta_file::fasta_file(std::string filename, open_type type, id_format format): bio_file(filename, type) {
    m_format = format;
}

std::unique_ptr <bio_entry> fasta_file::readLine() {
    char tchar('\0'), strand('\0');
    std::string sequence(""), chr(""), s_start(""), s_end("");
    long start(0), end(0);
    bool isID(true), continue_reading(true), isStart(false), isEnd(false), isStrand(false);
    if(getType() == write) {
        throw std::invalid_argument("Error : file opened in writing mode");
    }
    tchar = getChar(); // read first char to remove >
    while(remainToRead() && continue_reading) {
        tchar = getChar();
        if(tchar != '\n') {
            if(isID) { // if :\d-\d
                if(m_format == standard) {
                    chr += tchar;
                } else if(m_format == bedtools || m_format == bedtools_stranded) {
                    if(tchar == ':') {
                        isStart = true;
                    } else if(tchar == '-') {
                        isEnd = true;
                        isStart = false;
                    } else if(tchar == '(') {
                        isStrand = true;
                        isEnd = false;
                    } else if(tchar == ')') {
                        isStrand = false;
                        // nothing here : just there so as to pass the last parenthesis
                    } else if(isStart) {
                        s_start += tchar;
                    } else if(isEnd) {
                        s_end += tchar;
                    } else if(isStrand) {
                        strand = tchar;
                    } else { // none of these cases : its id
                        chr += tchar;
                    }
                }
            } else {
                sequence += toupper(tchar);
            }
        } else {
            if(isID) {
                isID = false;
            } else {
                if(peek() == '>') {
                    continue_reading = false;
                }
            }
        }
    }
    if(m_format == bedtools  || m_format == bedtools_stranded) {
        start = std::stol(s_start);
        end = std::stol(s_end);
    } else { // standard case
        start = 0;
        end = sequence.size();
        strand = '+';
    }
    if(m_format == bedtools) {
        strand = '+';
    }
    return std::make_unique <fasta_entry>(fasta_entry(chr, start, end, "none", strand, sequence));
    // return line;
}

std::string splitSequence (const std::string_view sequence, int limit = 80) {
    std::ostringstream new_seq;
    for(int i(0); i < (sequence.size() / 80 + 1); i++) {
        if(i % 100 == 0) {
            std::cout << i << "\r";
        }
        new_seq << sequence.substr(i * 80, 80);
        new_seq << '\n';
    }
    return new_seq.str();
}

std::string fasta_entry::getString() const {
    return std::string(">" + m_chr + "\n" + splitSequence(m_sequence));
}

fasta_entry* fasta_file::getEntry(int index) const {
    return dynamic_cast <fasta_entry*>(bio_file::getEntry(index));
}

char fasta_entry::getBase(int pos) const {
    try {
        return toupper(m_sequence.at(pos));
    } catch(std::out_of_range) {
        std::cout << "position " + std::to_string(pos) + " does not exist in " << m_chr << std::endl;
        std::cout << "max size : " << m_sequence.length() << std::endl;
        throw;
    }
}

void fasta_entry::mutate(char base, int pos) {
    try {
        m_sequence.at(pos);
        m_sequence[pos] = base;
    } catch(std::out_of_range) {
        std::cout << "Sequence size is " << m_sequence.size() << std::endl;
        std::cout << "Pos was " << pos << std::endl;
        std::cout << "Error : position specified is out of sequence bounds" << std::endl;
        throw;
    }
    
}

std::string reverseComp(std::string& sequence) {
    std::string sequence_final;
    for(const char& nuc: sequence) {
        switch(nuc) {
            case 'A':
                sequence_final = 'T' + sequence_final;
                break;
            case 'C':
                sequence_final = 'G' + sequence_final;
                break;
            case 'G':
                sequence_final = 'C' + sequence_final;
                break;
            case 'T':
                sequence_final = 'A' + sequence_final;
                break;
            default:
                sequence_final = nuc + sequence_final;
        }
    }
    return sequence_final;
}

std::string fasta_entry::subset(const bio_entry& mask) const {
    std::string return_string(m_sequence.substr(mask.getStart(), mask.getEnd() - mask.getStart()));
    if(mask.getStrand() == '-') {
        return_string = reverseComp(return_string);
    }
    return return_string;
}

std::string fasta_entry::subset(const int start, const int stop) const {
    return std::string(m_sequence.substr(start, stop - start));
}

void fasta_entry::apply_mask(std::vector <bio_entry*> mask) {
    std::string new_seq(m_sequence.size(), 'N');
    for(const bio_entry* value: mask) {
        new_seq.replace(value->getStart(), value->getSize(), subset(*value));
    }
    m_sequence = new_seq;
}

void fasta_file::apply_mask(const bio_file& file) {
    for(const std::string& chr: getChrs()) {
        std::cout << chr << "\r";
        fasta_entry* fasta_converted = dynamic_cast <fasta_entry*> (getEntriesByChr(chr)[0]);
        try {
            const std::vector <bio_entry*> mask(file.getEntriesByChr(chr));
            fasta_converted->apply_mask(mask);
        } catch(std::out_of_range) {
            fasta_converted->apply_mask(std::vector <bio_entry*>());
        }
    }
}

std::map<int, std::vector<std::shared_ptr <bio_entry>>> fasta_entry::matchPatterns(std::string pattern, bool complete_matchs) const {
    if(complete_matchs) {
        throw "not implemented";
    }
    // std::vector <std::unique_ptr <bio_entry>> results; 
    std::map<int, std::vector<std::shared_ptr <bio_entry>>> results; // 1 is int 2 is non-int
    if(pattern.size() == 2) {
        bool is_pattern(false);
        int start_non(0), start_true(0), phase_non(0), end_true(0);
        bool not_non(false), is_true(false);
        for(int i(0); i < m_sequence.size(); i++) {
            if(is_true) {
                if(m_sequence[i] != pattern[(i - start_true) % 2]) {
                    is_true = false;
                    end_true = i;
                    if((i - start_true) % 2 != 0) {
                        end_true --;
                    }
                    results[0].push_back(std::make_shared <bio_entry>(bio_entry(m_chr, start_true, end_true, '+', "CG")));
                }
                if(m_sequence[i] != 'N' & m_sequence[i] != pattern[(i - start_true) % 2]) {
                    not_non = false;
                    start_non = i;
                }
            } else if(not_non) {
                if(m_sequence[i] != 'N' & m_sequence[i] != pattern[(i - phase_non) % 2]) {
                    not_non = false;
                    start_non = i;
                }
            } else {
                if(m_sequence[i] == pattern[0] & m_sequence[i + 1] == pattern[1]) {
                    not_non = true;
                    is_true = true;
                    start_true = i; // set start
                }
                if((m_sequence[i] == pattern[0] | m_sequence[i] == 'N') & (m_sequence[i + 1] == pattern[1] | m_sequence[i + 1] == 'N')) {
                    not_non = true;
                    results[1].push_back(std::make_shared <bio_entry>(bio_entry(m_chr, start_non, i, '+', "nCG")));
                    phase_non = i;
                }
            }
        }
    } else {
        std::regex regex(pattern);
        std::sregex_iterator rit (m_sequence.begin(), m_sequence.end(), regex);
        std::sregex_iterator rend;

        std::string match(".");
        while(rit != rend) {
            const std::smatch smatch = *rit;
            int position = smatch.position();
            match = smatch.str();
            int length = position + match.size();
            results[0].push_back(std::make_shared <bio_entry>(bio_entry(m_chr, position, length, '+', match)));
            // bed_entry entry(m_header.getID(), position, length, match, 0, '.');
            // results.push_back(entry);
            rit ++;
        }
    }
    return results;
}

std::map<int, std::vector<std::shared_ptr <bio_entry>>> fasta_file::matchPatterns(std::string pattern, bool complete_matchs) const {
    std::map<int, std::vector<std::shared_ptr <bio_entry>>> results;
    std::map<int, std::vector<std::shared_ptr <bio_entry>>> tmp;
    for(const std::string& chr: getChrs()) {
        fasta_entry* fasta_converted = dynamic_cast <fasta_entry*> (getEntriesByChr(chr)[0]);
        std::cout << chr << "\r";
        tmp = fasta_converted -> matchPatterns(pattern, complete_matchs);
        results[0].insert(results[0].end(), tmp[0].begin(), tmp[0].end());
        results[1].insert(results[1].end(), tmp[1].begin(), tmp[1].end());
    }
    return results;
}

// std::unique_ptr<bio_entry> fasta_entry::intersect(const bio_entry& entry, bool stranded) {
//     intersect_results intersected_entry (bio_entry::intersect(entry, stranded)); // copy the value pointed
//     if(!(intersected_entry.resul == bio_entry())) {
//         return std::make_unique <fasta_entry> (fasta_entry(intersected_entry, subset(entry)));
//     } else {
//         return std::make_unique <fasta_entry> (fasta_entry());
//     }
// }