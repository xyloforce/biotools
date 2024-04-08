#pragma once
#include "bio_entry.hpp"
#include "bio_file.hpp"

class bed_entry: public bio_entry {
    public:
        bed_entry(std::string chr, long start, long end, std::string id, int score, char strand);
        bed_entry(bio_entry entry, int score);
        bed_entry();
        virtual ~bed_entry();

    // getters
        virtual std::string getString() const;
        int getScore() const;
    protected:
        int m_score;
};

class bed_file: public bio_file {
    public:
        bed_file(std::string filename, open_type type);

    // functions
        virtual std::unique_ptr <bio_entry> readLine();
        void appendVector(std::vector <std::shared_ptr <bio_entry>>& entries);
        void apply_intersect(bio_file& file, bool stranded = false, id_status status = both);

    private:
        // std::vector <std::unique_ptr<bed_entry>> m_content;
};

class AOE_entry: public bed_entry {
    public:
        AOE_entry(std::string chr, long start, long end, std::string id, int score, char strand, long zero);
        AOE_entry(bio_entry entry, int score, long zero);
        AOE_entry();

    // getters
        int getZero();
        int getRelativePos(int pos) const;
        virtual std::string getString() const;

    // functions
        // virtual std::unique_ptr<bio_entry> intersect(const bio_entry& entry, bool stranded);
        void appendVector(std::vector <std::shared_ptr <bio_entry>>& entries);

    private:
        long m_zero;
};

class AOE_file: public bio_file {
    public:
        AOE_file(std::string filename, open_type type);

    // functions
        virtual std::unique_ptr <bio_entry> readLine();
        void apply_intersect(bio_file& file, bool stranded = false, id_status status = both);
        // std::vector <AOE_entry> intersect(bio_file& file, bool stranded = false);

    private:
        // std::vector <std::unique_ptr<AOE_entry>> m_content;
};