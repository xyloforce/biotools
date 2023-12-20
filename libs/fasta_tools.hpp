#pragma once
#include "bio_file.hpp"
#include "bio_entry.hpp"

enum id_format {bedtools, bedtools_stranded, standard};

class fasta_entry: public bio_entry {
    public:
        fasta_entry(std::string chr, long start, long end, std::string id, char strand, std::string sequence);
        fasta_entry(bio_entry entry, std::string sequence);
        fasta_entry();

    // getters
        char getBase(int pos) const;
        virtual std::string getString() const;

    // functions
        void mutate(char base, int pos);
        std::string subset(const bio_entry& mask) const;
        std::string subset(const int start, const int stop) const;
        void apply_mask(std::vector <bio_entry*> mask);
        std::map<int, std::vector<std::shared_ptr <bio_entry>>> matchPatterns(std::string pattern, bool complete_matchs) const;
        // std::unique_ptr<bio_entry> intersect(const bio_entry& entry, bool stranded);

    private:
        std::string m_sequence;
};

class fasta_file: public bio_file {
    public:
        fasta_file(std::string filename, open_type type, id_format format);

    // functions
        virtual std::unique_ptr <bio_entry> readLine();
        void apply_mask(const bio_file& file);
        std::map<int, std::vector<std::shared_ptr <bio_entry>>> matchPatterns(std::string pattern, bool complete_matchs = false) const;
    
    // getters
        fasta_entry* getEntry(int index) const;

    private:
        id_format m_format;
};