#include "wiggle_tools.hpp"
#include "tools.hpp"
#include <cmath>
#include <iostream>

int main(int argc, char *argv[]) {
    std::cout << "Starting" << std::endl;
    wiggle_file test("test.wig", read);
    test.readWholeFile();
    std::cout << dynamic_cast <wiggle_entry*> (test.getEntry(1)) -> getString() << std::endl;
}