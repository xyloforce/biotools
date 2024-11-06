#include "bio_file.hpp"
#include <iostream>
#include <algorithm>
#include <sys/stat.h>

bio_file::bio_file(std::string filename, open_type type) {
    m_type = type;
    m_filename = filename;
    switch(type) {
        case read:
            openFile(filename, m_input_stream);
            break;
        case write:
            createFile(filename, m_output_stream);
            break;
    }
}

std::vector <std::string> bio_file::readBioLine(char delim) { // read one line and return an array of values
    char tchar('\0');
    std::vector <std::string> line(1, std::string());
    if(m_type == write) {
        std::cout << "Error : file opened in writing mode" << std::endl;
        exit(1);
    }
    while(tchar != '\n' && remainToRead()) {
        m_input_stream.get(tchar);
        if(tchar != '\n') {
            if(tchar != delim) {
                line[line.size() - 1] += tchar;
            } else {
                line.push_back(std::string());
            }
        }
    }
    return line;
}

void bio_file::readWholeFile() {
    std::cout << "Reading file " << m_filename << std::endl;
    bool warned = false;
    int i(0);
    while(remainToRead()) {
        try {
            if(i % 100 == 0) {
                std::cout << "Read " << i << " entries\r";
            }
            m_content.push_back(std::move(readLine()));
            m_indexes[m_content.at(i) -> getChr()].push_back(m_content.at(i).get()); // store a copy of the pointer for easy access
            i ++;
        } catch(const std::out_of_range& e) {
            if(!warned) {
                std::cout << "ignored line at position " << getLine() << std::endl;
                std::cout << "invalid line: " << e.what() << std::endl;
                warned = true;
            } else {
                std::cout << e.what() << std::endl;
                std::cout << "error happened at line " << getLine() << std::endl;
                std::cout << "more than one empty line: format is probably incorrect" << std::endl;
                throw;
            }
        }
    }
}

void bio_file::writeBioLine(const bio_entry& entry) {
    if(m_type != write) {
        std::cout << "Error : file opened in reading mode" << std::endl;
        exit(1);
    }
    m_output_stream << entry.getString() << "\n";
}

void bio_file::openFile(std::string filename, std::ifstream& stream) {
    struct stat buffer;
    if(stat(filename.c_str(), &buffer) != 0) {
        std::cout << "file " << filename << " does not exist" << std::endl;
        exit(1);
    }
    stream.open(filename);
} // open the file in reading mode

void bio_file::createFile(std::string filename, std::ofstream& stream) {
    stream.open(filename);
} // open the file in writing mode

void bio_file::writeToFile() {
    for(const auto& entry: m_content) {
        writeBioLine(*entry);
    }
}

bool bio_file::remainToRead() const {
    return !m_input_stream.eof();
}

std::vector <intersect_results> bio_file::intersect(const bio_file& file, bool stranded, id_status status) {
    std::vector <intersect_results> results;
    for(const auto& pair: m_indexes) { // pair chr : pointers
        try {
            const std::vector <bio_entry*> vals = file.getEntriesByChr(pair.first);
            int A(0), B(vals.size() - 1);
            bool found(false);
            int index2(0);
            for(const auto& entry: pair.second) {
                found = false;
                A = 0;
                B = vals.size() - 1;
                while(A <= B && !found) {
                    index2 = int((A+B)/2);
                    intersect_results intersected_entry(entry -> intersect(vals[index2], stranded, status));
                    if(intersected_entry.result == bio_entry()) {
                        if((*entry) < *(vals[index2])) {
                            B = index2 - 1;
                        } else {
                            A = index2 + 1;
                        }
                    } else { // we got a hit, is it the only hit ?
                        results.push_back(intersected_entry);
                        for(int i(index2 - 1); i > 0; i --) {
                            intersect_results intersected_entry(entry -> intersect(vals[i], stranded, status));
                            if(intersected_entry.result == bio_entry()) {
                                break;
                            } else {
                                results.push_back(intersected_entry);
                            }
                        }
                        for(int i(index2 + 1); i < vals.size(); i ++) {
                            intersect_results intersected_entry(entry -> intersect(vals[i], stranded, status));
                            if(intersected_entry.result == bio_entry()) {
                                break;
                            } else {
                                results.push_back(intersected_entry);
                            }
                        }
                        found = true;
                    }
                }
            }
        } catch(std::out_of_range) {
            std::cout << "no entries for this chr, skipping" << std::endl;
        }
    }
    return results;
}

void bio_file::apply_intersect(bio_file& file, bool stranded, id_status status) {
    std::vector <intersect_results> results(intersect(file, stranded, status));
    clear(); // ATTENTION ICI ON VIDE LE TABLEAU DONC SOURCE NE POINTE PLUS VERS RIEN
    for(const auto& entry: results) {
        appendEntry(std::make_unique <bio_entry> (entry.result));
    }
}

bool true_sort(const bio_entry* A, const bio_entry* B) {
    return *A < *B;
}

std::vector <bio_entry*> bio_file::getEntriesByChr(const std::string chr) const {
    std::vector <bio_entry*> results;
    try {
        results = m_indexes.at(chr);
        sort(results.begin(), results.end(), true_sort);
    } catch(std::out_of_range) {
        std::cout << "No entry matching " + chr << std::endl;
        throw;
    }
    return results;
}

open_type bio_file::getType() const {
    return m_type;
}

char bio_file::getChar() {
    char returned('\0');
    m_input_stream.get(returned);
    return returned;
}

bio_entry* bio_file::getEntry(int index) const {
    return m_content.at(index).get();
}

int bio_file::getSize() const {
    return m_content.size();
}

int bio_file::getLine() {
    if(m_type == read) {
        int current_pos(m_input_stream.tellg()), number_lines(0);
        if(m_input_stream.eof()) {
            return -1;
        } else if(m_input_stream.fail()) {
            std::cout << "issues with the stream, investigate" << std::endl;
            return -1;
        } else {
            char tchar('\0');
            m_input_stream.seekg(0, std::ios::beg); // go back to the begining
            while(m_input_stream.tellg() < current_pos) {
                m_input_stream.get(tchar);
                if(tchar == '\n') {
                    number_lines ++;
                }
            }
            return number_lines;
        }
    } else {
        throw std::invalid_argument("not implemented for file writing");
    }
}

void bio_file::writeString(std::string value) {
    if(m_type == write) {
        m_output_stream << value;
    } else {
        throw std::domain_error("tried a write operation in a read-only file");
    }
}

void bio_file::appendEntry(std::unique_ptr <bio_entry> entry) {
    m_content.push_back(std::move(entry));
    m_indexes[m_content.back() -> getChr()].push_back(m_content.back().get());
}

void bio_file::typeToWrite(const std::string filename) {
    m_input_stream.close();
    createFile(filename, m_output_stream);
    m_type = write;
    m_filename = filename;
}

std::vector <std::string> bio_file::getChrs() const {
    std::vector <std::string> results;
    for(const auto& entry: m_indexes) {
        results.push_back(entry.first);
    }
    return results;
}

char bio_file::peek() {
    return m_input_stream.peek();
}

void bio_file::clear() {
    m_content.clear();
    m_indexes.clear();
}

std::vector <bio_entry*> bio_file::getEntries() const {
    std::vector <bio_entry*> results;
    for(const auto& entry: m_content) {
        results.push_back(entry.get());
    }
    return results;
}

void bio_file::eraseAndLoadBlock(int amount) {
    int amount_loaded(0);
    clear();
    while(remainToRead() && amount > amount_loaded) {
        try  {
            appendEntry(readLine());
        } catch(const std::out_of_range& e) {
            std::cout << "ignored line at position " << getLine() << std::endl;
            std::cout << "invalid line: " << e.what() << std::endl;
        }
        amount_loaded ++;
        if(amount_loaded % 100 == 0) {
            std::cout << "Loaded " << amount_loaded << " entries\r";
        }
    }
}