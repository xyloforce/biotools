#include "tools.hpp"
#include <algorithm>
#include <iostream>

std::map <char, std::string> getArgs(std::vector <std::string> args) {
    args.erase(args.begin(), args.begin() + 1); // delete first element (program name)
    std::map <char, std::string> parsed;
    char last_arg('\0');
    for(const std::string& arg: args) {
        if(arg[0] == '+') { // its an argument
            parsed[arg[1]] = "";
            last_arg = arg[1];
        } else {
            parsed[last_arg] = arg;
        }
    }
    return parsed;
}