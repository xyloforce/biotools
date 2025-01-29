#include "vcf_tools.hpp"
#include <iostream>

vcf_file::vcf_file(std::string filename, open_type type): bio_file(filename, type) {
    if(type == write) {
        writeString("##fileformat=VCFv4.2\n");
        writeString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    } else if(type == read) {
        readBioLine('\t');
        readBioLine('\t');
    } else if(type == none) {
        // nothing bc the point is that it's unset
    }
}

vcf_file::vcf_file(): bio_file() {

}

void vcf_file::typeToWrite(const std::string filename) {
    bio_file::typeToWrite(filename);
    writeString("##fileformat=VCFv4.2\n");
    writeString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
}

vcf_entry::vcf_entry(std::string chr, long start, long end, std::string id, std::string ref, std::vector <std::string> alt, int qual, std::string filter, std::map <std::string, std::string> infos): bio_entry(chr, start, end, '+', id) {
    if(ref.size() > 0 && alt.size() > 0) {
        m_ref = ref;
        m_alt = alt;
        m_qual = qual;
        m_filter = filter;
        m_infos = infos;
    } else {
        std::invalid_argument("no value for either alt or ref, probable parsing issue");
    }
}

vcf_entry::vcf_entry(bio_entry entry, std::string ref, std::vector <std::string> alt, int qual, std::string filter, std::map <std::string, std::string> infos): bio_entry(entry), m_ref(ref), m_alt(alt), m_qual(qual), m_filter(filter), m_infos(infos) {

}

vcf_entry::vcf_entry(): bio_entry() {
    m_ref = "";
    m_alt = std::vector <std::string>();
    m_qual = 0;
    m_filter = "";
    m_infos = std::map <std::string, std::string>();
}

std::string vcf_entry::getString() const {
    std::string s_alt, s_info;
    for(const std::string& value: m_alt) {
        s_alt += value + ',';
    }
    if(s_alt.size() > 1) {
        s_alt.pop_back();
    } // remove the last comma
    for(const auto& pair: m_infos) {
        s_info += pair.first + '=' + pair.second + ';';
    }
    if(s_info.size() > 1) {
        s_info.pop_back();
    }
    if(s_alt == "") {
        s_alt = ".";
    }
    if(s_info == "") {
        s_info = ".";
    }
    return m_chr + "\t" + std::to_string(m_start + 1) + "\t" + m_id + "\t" + m_ref + "\t" + s_alt + "\t" + std::to_string(m_qual) + "\t" + m_filter + "\t" + s_info;
    // return m_chr + "\t" + std::to_string(m_start) + "\t" + m_id + "\t" + m_ref + "\t" + s_alt + "\t" + std::to_string(m_qual) + "\t" + m_filter + "\t" + s_info;
}

std::unique_ptr <bio_entry> vcf_file::readLine() {
    std::string chr, id, ref, filter, t_key, t_value;
    long start, end;
    std::vector <std::string> alt(1, std::string());
    std::map <std::string, std::string> infos;
    char tchar('\0');
    int qual(0);
    bool isKey(true);
    std::vector <std::string> line = readBioLine('\t');
    if(line.at(0)[0] != '#') {
        if(line.size() >= 8) {
            chr = line.at(0);
            try {
                start = std::stol(line.at(1)) - 1;
            } catch(std::invalid_argument) {
                std::cout << "Invalid conversion : " << line.at(1) << "\r";
                throw;
            }
            end = start + 1;
            id = line.at(2);
            ref = line.at(3);
            for(int i(0); i < line.at(4).size(); i++) {
                tchar = line.at(4)[i];
                if(tchar != ';' && tchar != ' ') {
                    alt[alt.size() - 1] += tchar;
                }
            }
            if(line.at(5) != ".") {
                qual = std::stoi(line.at(5));       
            } else {
                qual = 0;
            }
            filter = line.at(6);
            for(int i(0); i < line.at(7).size(); i++) {
                tchar = line.at(7)[i];
                if(tchar != ' ') { // skip whitespaces
                    if(tchar != ';') { // if it is, need to reset tkey/tval
                        if(tchar != '=') { // change from key to val
                            if(isKey) {
                                t_key += tchar;
                            } else {
                                t_value += tchar;
                            }
                        } else {
                            isKey = false;
                        }
                    } else {
                        infos[t_key] = t_value;
                        t_key = "";
                        t_value = "";
                        isKey = true;
                    }
                }
            }
            return std::make_unique <vcf_entry> (vcf_entry(chr, start, end, id, ref, alt, qual, filter, infos));
        } else {
            throw(std::out_of_range("line doesn't have the right number of elements : " + std::to_string(line.size())));
        }
    } else {
        throw(std::invalid_argument("skipped comment line"));
    }
}

std::string vcf_entry::getRef() const {
    return m_ref;
}

std::string vcf_entry::getAlt(int number) const {
    try {
        return m_alt.at(number);
    } catch(std::out_of_range) {
        std::cout << "Asked for alt value number " << number << " where only " << m_alt.size() << " are available\n";
        std::cout << "No base available at this position" << std::endl;
        throw;
    }
}

// std::unique_ptr<bio_entry> vcf_entry::intersect(const bio_entry& entry, bool stranded) {
//     bio_entry intersected_entry (*bio_entry::intersect(entry, stranded)); // copy the value pointed
//     if(!(intersected_entry == bio_entry())) {
//         return std::make_unique <vcf_entry> (vcf_entry(intersected_entry, m_ref, m_alt, m_qual, m_filter, m_infos));
//     } else {
//         return std::make_unique <vcf_entry> (vcf_entry());
//     }
// }

void vcf_file::eraseAndLoad() {
    clear();
    try {
        appendEntry(std::move(readLine()));
    } catch(std::out_of_range) {
        std::cout << "is an empty line, probably eof" << std::endl;
    }
    
}
