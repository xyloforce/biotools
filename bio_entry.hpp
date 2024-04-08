#pragma once
#include <string>
#include <memory>

struct intersect_results;
enum id_status {source, hit, both};

class bio_entry {
    public:
        bio_entry(std::string chr, long start, long end, char strand, std::string id);
        bio_entry();
        virtual ~bio_entry();

        // getters
        std::string getChr() const;
        long getStart() const;
        long getEnd() const;
        char getStrand() const;
        std::string getID() const;
        int getSize() const;

        // setters
        void setID(std::string id);
        void setStart(int start);
        void setEnd(int end);
        
        // operators
        bool operator == (const bio_entry& entry) const;
        bool operator < (const bio_entry& entry) const;

        // functions
        virtual std::string getString() const;
        virtual intersect_results intersect(const bio_entry *entry, bool stranded, id_status status = both);

    protected:
        std::string m_chr;
        long m_start;
        long m_end;
        char m_strand;
        std::string m_id;
};

struct intersect_results {
    bio_entry result;
    bio_entry* source;
    const bio_entry* hit;
};