#include "bed_tools.hpp"
#include <iostream>

bed_entry::bed_entry(std::string chr, long start, long end, std::string id, int score, char strand): bio_entry(chr, start, end, strand, id) {
    m_score = score;
}

bed_entry::bed_entry(bio_entry entry, int score): bio_entry(entry) {
    m_score = score;
}

bed_entry::bed_entry(): bio_entry() {
    m_score = 0;
}

bed_entry::~bed_entry() {
    
}

std::string bed_entry::getString() const {
    return m_chr + "\t" + std::to_string(m_start) + "\t" + std::to_string(m_end) + "\t" + m_id + "\t" + std::to_string(m_score) + "\t" + std::string(1, m_strand);
}

bed_file::bed_file(std::string filename, open_type type) : bio_file(filename, type) {

}

std::unique_ptr <bio_entry> bed_file::readLine() {
    bed_entry entry;
    if(remainToRead()) {
        std::vector <std::string> values = readBioLine('\t');
        try {
            std::string chr = values.at(0);
            long start = std::stol(values.at(1));
            long end = std::stol(values.at(2)); // minimally defined bed entry
            std::string name;
            int score(0);
            char strand;
            try {
                name = values.at(3);
            } catch(const std::out_of_range) {
                name = ".";
            }
            try {
                score = std::stoi(values.at(4));
            } catch(std::logic_error) { // catch both score undefined (= . so no conversion) and no score col
                score = 0;
            }
            try {
                strand = values.at(5)[0];
            } catch(const std::out_of_range) {
                strand = '+';
            } // these cols are optionnal
            entry = bed_entry(chr, start, end, name, score, strand);
        }
        catch(const std::out_of_range)
        {
            std::cout << "not enough elements in line" << std::endl;
            throw;
        }
    }
    return std::make_unique <bed_entry>(entry);
}

void bed_file::apply_intersect(bio_file& file, bool stranded, id_status status) {
    std::vector <intersect_results> results(intersect(file, stranded, status));
    clear();
    for(const auto& entry: results) {
        appendEntry(std::make_unique <bed_entry> (bed_entry(entry.result, 1000)));
    }
}

AOE_entry::AOE_entry(std::string chr, long start, long end, std::string id, int score, char strand, long zero): bed_entry(chr, start, end, id, score, strand) {
    m_zero = zero;
}

AOE_entry::AOE_entry(bio_entry entry, int score, long zero): bed_entry(entry, score) {
    m_zero = zero;
}

AOE_entry::AOE_entry(): bed_entry() {
    m_zero = 0;
}

int AOE_entry::getZero() {
    return m_zero;
}

std::string AOE_entry::getString() const {
    return m_chr + "\t" + std::to_string(m_start) + "\t" + std::to_string(m_end) + "\t" + m_id + "\t" + std::to_string(m_score) + "\t" + std::string(1, m_strand) + "\t" + std::to_string(m_zero);
}

AOE_file::AOE_file(std::string filename, open_type type): bio_file(filename, type) {

}

std::unique_ptr <bio_entry> AOE_file::readLine() {
    AOE_entry entry;
    if(remainToRead()) {
        std::vector <std::string> values = readBioLine('\t');
        try {
            char strand = '\0';
            if(values.at(5)[0] == 'L') {
                strand = '+';
            } else if(values.at(5)[0] == 'R') {
                strand = '-';
            } else if(values.at(5)[0] == '+' || values.at(5)[0] == '-') {
                strand = values.at(5)[0];
            } else {
                throw std::domain_error("strand incorrectly defined : not a correct value : " + std::string(1, values.at(5)[0]));
            } // correct strand info for AOE files
            int score(0);
            try {
                score = std::stoi(values.at(4));
            } catch(std::invalid_argument) {
                score = 0;
            }
            entry = AOE_entry(values.at(0), std::stol(values.at(1)), std::stol(values.at(2)), values.at(3), score, strand, std::stol(values.at(6)));
        }
        catch(const std::out_of_range)
        {
            std::cout << "not enough elements in line" << std::endl;
            throw;
        }
    }
    return std::make_unique <AOE_entry>(entry);
}

int AOE_entry::getRelativePos(int pos) const {
    int returned(pos - m_zero);
    if(m_strand == '-') {
        returned = returned * -1;
    }
    return returned;
}

// std::unique_ptr<bio_entry> AOE_entry::intersect(const bio_entry& entry, bool stranded) {
//     bio_entry intersected_entry (*bio_entry::intersect(entry, stranded)); // copy the value pointed
//     if(!(intersected_entry == bio_entry())) {
//         return std::make_unique <AOE_entry> (AOE_entry(intersected_entry, 1000, m_zero));
//     } else {
//         return std::make_unique <AOE_entry> (AOE_entry());
//     }
// }

void bed_file::appendVector(std::vector <std::shared_ptr <bio_entry>>& entries) {
    for(const auto& entry: entries) {
        bed_entry b_entry((*entry), 1000);
        appendEntry(std::make_unique <bed_entry>(b_entry));
    }
}

void AOE_file::apply_intersect(bio_file& file, bool stranded, id_status status) {
    std::vector <intersect_results> results(intersect(file, stranded, status));
    std::vector <std::unique_ptr <bio_entry>> entries;
    for(const auto& entry: results) {
        entries.push_back(std::make_unique <AOE_entry> (AOE_entry(entry.result, dynamic_cast <AOE_entry*>(entry.source) -> getScore(), dynamic_cast <AOE_entry*>(entry.source) -> getZero())));
    }
    clear();
    for(auto &entry: entries) {
        appendEntry(std::move(entry));
    }
}

int bed_entry::getScore() const {
    return m_score;
}
