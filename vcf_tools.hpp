#pragma once
#include "bio_file.hpp"
#include "bio_entry.hpp"

class vcf_entry: public bio_entry {
    public:
        vcf_entry(std::string chr, long start, long end, std::string id, std::string ref, std::vector <std::string> alt, int qual, std::string filter, std::map <std::string, std::string> infos);
        vcf_entry(bio_entry entry, std::string ref, std::vector <std::string> alt, int qual, std::string filter, std::map <std::string, std::string> infos);
        vcf_entry();

    // getters
        virtual std::string getString() const;
        std::string getRef() const;
        std::string getAlt(int number = 0) const;

    // functions
        // virtual std::unique_ptr<bio_entry> intersect(const bio_entry& entry, bool stranded);

    private:
        std::string m_ref;
        std::vector <std::string> m_alt;
        int m_qual;
        std::string m_filter;
        std::map <std::string, std::string> m_infos;
};

class vcf_file: public bio_file {
    public:
        vcf_file(std::string filename, open_type type);
        vcf_file();
    // functions
        virtual std::unique_ptr <bio_entry> readLine();
        void typeToWrite(std::string filename);
        void eraseAndLoad();
};