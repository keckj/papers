
#include "headers.hpp"
#include <iostream>
#include "rand.hpp"

#include "interval.hpp"
#include "integerGrid.hpp"
#include "binaryTreeNode.hpp"
#include "sample.hpp"
#include "functionSample.hpp"
#include "point.hpp"
#include "interval.hpp"
#include "gnuplot.hpp"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
#include <conio.h>   //for getch(), needed in wait_for_key()
#include <windows.h> //for Sleep()
void sleep(int i) { Sleep(i*1000); }
#endif

void wait_for_key();

void header();
void footer();

float F(float t) {
    return sin(10*t);
}

int main(int argc, char **argv) {

    header();

    Random::init();

    Globals::init();
    Globals::check();
    
    Gnuplot gp("tee plot.gp | gnuplot -persist");

    Interval<float> unitInterval(0.0f,1.0f);
    Interval<float> halfInterval(0.25f,0.75f);

    BinaryTreeNode<float> *tree = new BinaryTreeNode<float>(1,1,0.0,1.0);

    FunctionSample<1000u, float> sample(halfInterval, F);
    for (unsigned int i = 0; i < 1000; i++) {
        tree->insert(sample[i],2);
    }
    //sample.plot(gp);

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

void wait_for_key ()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
    std::cout << std::endl << "Press any key to continue..." << std::endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    std::cout << std::endl << "Press ENTER to continue..." << std::endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}
