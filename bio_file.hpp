#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <map>
#include "bio_entry.hpp"

enum open_type {read, write, none};

class bio_file {
    public:
        bio_file(std::string filename, open_type type);
        bio_file();
    // getters
        std::vector <bio_entry*> getEntriesByChr(const std::string chr) const;
        std::vector <bio_entry*> getEntries() const; 
        bio_entry* getEntry(int index) const;
        open_type getType() const;
        std::vector <std::string> getChrs() const;
        char getChar();
        int getSize() const; // return size in number of entries
        int getLine();
        char peek();

    // functions
        std::vector <std::string> readBioLine(char delim); // return a vector containing the components
        virtual std::unique_ptr <bio_entry> readLine() = 0; // return a pointer towards a bioentry 
        void writeBioLine(const bio_entry& entry); // write the string representation of a bioentry
        void openFile(std::string filename, std::ifstream& stream); // open the file in reading mode
        void createFile(std::string filename, std::ofstream& stream); // open the file in writing mode
        virtual void readWholeFile();
        virtual void writeToFile(); // write a full array of bio_entries
        bool remainToRead() const;
        std::vector <intersect_results> intersect(const bio_file& file, bool stranded = false, id_status status = both);
        virtual void apply_intersect(bio_file& file, bool stranded = false, id_status status = both);
        void writeString(std::string value);
        virtual void appendEntry(std::unique_ptr <bio_entry> entry);
        void typeToWrite(const std::string filename = "");
        void typeToRead(const std::string filename);
        void clear();
        void eraseAndLoadBlock(int amount = 10000);

    private:
        std::vector <std::unique_ptr<bio_entry>> m_content;
        std::map <std::string, std::vector <bio_entry*>> m_indexes;
        std::ifstream m_input_stream;
        std::ofstream m_output_stream;
        open_type m_type;
        std::string m_filename;
};
