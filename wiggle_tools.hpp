#pragma once
#include "bio_file.hpp"
#include "bio_entry.hpp"

struct valid_double;

class wiggle_entry: public bio_entry {
    public:
        wiggle_entry(std::string chr, long start, long end, std::string id, char strand, std::vector <valid_double> values, std::vector<std::map <std::string, std::string>> track_parameters);

        // getters
            std::string getString() const;
            std::map <std::string, std::string> getParameters(int line) const;
            std::vector <valid_double> getSubset(const bio_entry* entry) const;

private:
    std::vector <valid_double> m_values;
    std::vector<std::map <std::string, std::string>> m_track_parameters;
};

class wiggle_file: public bio_file {
    public:
        wiggle_file(std::string filename, open_type type);

    // functions
        virtual std::unique_ptr <bio_entry> readLine();
};

struct valid_double {
    valid_double();
    valid_double(double i_value, bool i_is_valid);
    double value;
    bool is_valid;
};