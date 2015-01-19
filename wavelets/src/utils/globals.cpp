
#include "headers.hpp"

#include "globals.hpp"
#include "utils.hpp"

#include <iostream>

void Globals::init() {
    std::cout << "[Globals] Init" << std::endl;
}

void Globals::check() {
}

void Globals::print(std::ostream &out) { 
    std::cout << "[Globals]" << std::endl;
    std::cout << "\t" << std::endl;
}

