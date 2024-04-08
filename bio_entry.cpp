#include "bio_entry.hpp"
#include <iostream>

bio_entry::bio_entry(std::string chr, long start, long end, char strand, std::string id): m_chr(chr), m_start(start), m_end(end), m_strand(strand), m_id(id) {
    if(start > end) {
        std::cout << "parameters were : " << chr << " " << start << " " << end << " " << strand << " " << id << "\n";
        throw std::invalid_argument("start > end, incorrect definition");
    }
}

bio_entry::bio_entry(): m_chr("undefined"), m_start(0), m_end(0), m_strand('\0'), m_id("none") {
    
}

bio_entry::~bio_entry() {
    
}

std::string bio_entry::getChr() const {
    return m_chr;
}

long bio_entry::getStart() const {
    return m_start;
}

long bio_entry::getEnd() const {
    return m_end;
}

char bio_entry::getStrand() const {
    return m_strand;
}

std::string bio_entry::getID() const {
    return m_id;
}

std::string bio_entry::getString() const {
    return m_chr + "_" + std::to_string(m_start) + "_" + std::to_string(m_end) + "_" + std::string(1, m_strand);
}

bool bio_entry::operator<(const bio_entry& entry) const {
    if(entry.getChr() == m_chr) {
        if(m_start == entry.getStart()) {
            if(m_end == entry.getEnd()) {
                return m_id < entry.getID();
            } else {
                return m_end < entry.getEnd();
            }
        } else {
            return m_start < entry.getStart();
        }
        return m_start < entry.getStart();
    } else {
        return m_chr < entry.getChr();
    }
}

bool bio_entry::operator==(const bio_entry& entry) const {
    return m_chr == entry.getChr() &&
           m_start == entry.getStart() &&
           m_end == entry.getEnd() &&
           m_strand == entry.getStrand() &&
           m_id == entry.getID();
}

intersect_results bio_entry::intersect(const bio_entry *entry, bool stranded) {
    bio_entry result;
    if(entry -> getChr() == m_chr) {
        if(entry -> getStrand() == m_strand || !stranded) {
            if(entry -> getStart() >= m_start && entry -> getEnd() <= m_end) {
                result = bio_entry(m_chr, entry -> getStart(), entry -> getEnd(), m_strand, m_id + "_" + entry -> getID());
            } else if(entry -> getStart() >= m_start && entry -> getStart() < m_end && entry -> getEnd() > m_end) {
                result = bio_entry(m_chr, entry -> getStart(), m_end, m_strand, m_id + "_" + entry -> getID());
            } else if(entry -> getStart() < m_start && entry -> getEnd() > m_start && entry -> getEnd() <= m_end) {
                result = bio_entry(m_chr, m_start, entry -> getEnd(), m_strand, m_id + "_" + entry -> getID());
            } else if(entry -> getStart() <= m_start && entry -> getEnd() > m_start && entry -> getStart() < m_end && entry -> getEnd() >= m_end) {
                result = bio_entry(m_chr, m_start, m_end, m_strand, m_id + "_" + entry -> getID());
            }
        }
    }
    intersect_results results;
    results.result = result;
    results.source = this;
    results.hit = entry; 
    return results;
}

int bio_entry::getSize() const {
    return m_end - m_start;
}