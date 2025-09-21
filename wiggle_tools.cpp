#include "wiggle_tools.hpp"
#include <iostream>
#include <sstream>

wiggle_entry::wiggle_entry(std::string chr, long start, long end, std::string id, char strand, std::vector <valid_double> values, std::vector<std::map <std::string, std::string>> track_parameters): bio_entry(chr, start, end, strand, id), m_values(values), m_track_parameters(track_parameters) {

}

wiggle_file::wiggle_file(std::string filename, open_type type): bio_file(filename, type) {

}

std::unique_ptr <bio_entry> wiggle_file::readLine() {
    char tchar;
    int past_pos(0), step_entry(0);
    std::string key, value, s_pos, s_value, chr, id;
    bool is_track(true), is_key(true), continue_reading(true), is_def(false), is_fixed(false), is_inside_quotes(false);
    std::vector<std::map <std::string, std::string>> infos(2, std::map <std::string, std::string>());
    std::vector <valid_double> values;
    long start(0);

    while(remainToRead() && continue_reading) {
        tchar = getChar();
        if(is_track) { // first there is a track IF THERE IS NO TRACK ADD IT
            if((tchar == ' ' && !is_inside_quotes) || tchar == '\n') {
                is_key = true;
                if(value != "") {
                    infos[0][key] = value;
                    key = "";
                    value = "";
                }
                if(tchar == '\n') {
                    is_track = false;
                    is_def = true;
                }
            } else if(tchar == '"') {
                is_inside_quotes = !is_inside_quotes;
            } else if(tchar == '=') {
                is_key = false;
            } else if(is_key) {
                key += tchar;
            } else {
                value += tchar;
            }
        } else if(is_def) { // its followed by a def line
            if (tchar == ' ' || tchar == '\t' || tchar == '\n') {
                if(key == "fixedStep" || key == "variableStep") {
                    is_fixed = (key == "fixedStep");
                    is_key = true;
                    key = "";
                } else {
                    is_key = true;
                    infos[1][key] = value;
                    key = "";
                    value = "";
                }
                if(tchar == '\n') {
                    is_def = false;
                    if(is_fixed) {
                        try {
                            start = std::stol(infos[1].at("start")) - 1; // wigs are 1-based
                        } catch(std::out_of_range) {
                            std::cout << "Error : start was not parsed" << std::endl;
                            throw;
                        }
                        try {
                            step_entry = std::stol(infos[1].at("step"));
                        } catch(std::out_of_range) {
                            step_entry = 1;
                        }
                    } else {
                        is_key = true;
                    }
                }
            } else if(tchar == '=') {
                is_key = false;
            } else if(is_key) {
                key += tchar;
            } else {
                value += tchar;
            }
        } else { // not a def or a track : either it's a numeric value or the start of a new int
            if(tchar == '\n') {
                try {
                    if(peek() == 't') { // if a step has no track header FUCKIN FIX IT
                        continue_reading = false;
                    }
                    /*if(std::stod(s_value) > 1) {
                        std::cout << "error : " << s_value << std::endl;
                        throw(666);
                    } why the hell did i wrote that here ???? */
                    if(!is_fixed) {
                        if(past_pos != 0) {
                            for(int i(0); i < (std::stoi(s_pos) - past_pos - 1); i++) { // same shit than below
                                values.push_back(valid_double());
                            }
                        } else {
                            start = std::stol(s_pos) - 1; // wig are 1-based
                        }
                        past_pos = std::stol(s_pos);
                        is_key = true;
                    } else {
                        if(values.size() != 0) {
                            for(int i(0); i < (step_entry - 1); i++) { // span of 100 : need to add values between 101 and 199 so 198 values
                                values.push_back(valid_double(std::stod(s_value), true));
                            }
                        }
                    }
                    values.push_back(valid_double(std::stod(s_value), true));
                } catch(std::invalid_argument) {
                    std::cout << "Conversion error : tried to convert " << s_value << " or " << s_pos << " to number." << std::endl;
                    throw;
                }
                s_value = "";
                s_pos = "";
            } if(!is_fixed) {
                if(tchar == ' ' || tchar == '\t') {
                    is_key = false;
                } else if(is_key) {
                    s_pos += tchar;
                } else {
                    s_value += tchar;
                }
            } else {
                s_value += tchar;
            }
        }
    } // file is track then fixed or variable then data
    try {
        chr = infos[1].at("chrom");
    } catch(std::out_of_range) {
        std::cout << "Error : chr wasnt parsed" << std::endl;
        throw;
    }
    return std::make_unique <wiggle_entry> (wiggle_entry(chr, start, start + values.size(), id, '+', values, infos));
}

std::string wiggle_entry::getString() const {
    std::ostringstream stream;
    bool is_fixed(true);
    int pos(getStart() + 1); // wig are 1-based but loaded as 0-based for easy intersection
    
    stream << "track ";
    for(const auto& pair: getParameters(0)) {
        stream << pair.first << "=\"" << pair.second << "\" ";
    }
    stream.seekp(long(stream.tellp()) - 1);
    stream << '\n';
    
    for(const auto& value: m_values) {
        if(!value.is_valid) {
            is_fixed = false;
            break;
        }
    }
    if(is_fixed) {
        stream << "fixedStep ";
    } else {
        stream << "variableStep ";
    }
    
    for(const auto& pair: getParameters(1)) {
        stream << pair.first << "=\"" << pair.second << "\" ";
    }
    stream.seekp(long(stream.tellp()) - 1);
    stream << "\n";
    for(int i(0); i < m_values.size(); i++) {
        if(is_fixed) {
            if(i % stoi(getParameters(1)["step"]) == 0) {
                stream << m_values[i].value << '\n';
            }
        } else {
            if(m_values[i].is_valid) {
                stream << pos << " " << m_values[i].value << "\n";
            }
            pos ++;
        }
    }
    return stream.str();
}

std::map <std::string, std::string> wiggle_entry::getParameters(int line) const {
    return m_track_parameters.at(line);
}

std::vector <valid_double> wiggle_entry::getSubset(const bio_entry* entry) const {
    std::vector <valid_double> results;
    if(getStart() - entry -> getStart() < 0 & entry -> getEnd() - getStart() - 1 >= m_values.size()) {
        std::cout << "values of entry are : " << entry -> getStart() << " - " << entry -> getEnd() << std::endl;
        std::cout << "values of current are : " << getStart() << " - " << getEnd() << std::endl;
        throw std::out_of_range("invalid range");
    }
    int range_start(getStart() - entry -> getStart());
    for(int i(entry -> getStart() - getStart()); i < (entry -> getEnd() - getStart() - 1); i ++) {
        results.push_back(m_values.at(i));
    }
    // return std::vector <valid_double>(&m_values[getStart() - entry -> getStart()], &m_values[entry -> getEnd() - getStart() - 1]);
    return results;
}

valid_double::valid_double(): value(0.0), is_valid(false) {}

valid_double::valid_double(double i_value, bool i_is_valid): value(i_value), is_valid(i_is_valid) {}