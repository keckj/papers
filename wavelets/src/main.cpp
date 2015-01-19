
#include "headers.hpp"
#include <iostream>
#include "rand.hpp"

#include "interval.hpp"
#include "integerGrid.hpp"
#include "binaryTreeNode.hpp"

void header();
void footer();

int main(int argc, char **argv) {

    header();
    
    Random::init();

    Globals::init();
    Globals::check();
    
    Interval<float> unitInterval(0.0f,1.0f);
    
    BinaryTreeNode<float> *node = new BinaryTreeNode<float>(0,0,0.0,1.0);
    std::cout << *node << std::endl;
    
    footer();

    return EXIT_SUCCESS;
}

void header() {
    std::cout << std::endl;
    std::cout << "*********************************************" << std::endl;
    std::cout << "**     Wavelet implementation v. " STRINGIFY_MACRO(WAVELETS_VERSION) "     **" << std::endl;
    std::cout << "*********************************************" << std::endl;
    std::cout << std::endl;
    std::cout << "Author:   Keck Jean-Baptiste -- Ensimag - MSIAM 2014-2015" << std::endl;
    std::cout << std::endl;
    std::cout << "Overview: This is a partial implementation of Christophe P. Bernard thesis on interpolating Delaurier and Dubuc wavelets in the context of wavelet paper review." << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Program is running in " STRINGIFY_MACRO(COMPILE_MODE) " mode !" << std::endl;
}

void footer() {
    std::cout << std::endl;
    std::cout << "All done, exiting !" << std::endl;
}
