
#include "headers.hpp"

#include <iostream>
#include <Eigen/Dense>

#include "rand.hpp"
#include "interval.hpp"
#include "integerGrid.hpp"
#include "binaryTreeNode.hpp"
#include "sample.hpp"
#include "functionSample.hpp"
#include "point.hpp"
#include "interval.hpp"
#include "gnuplot.hpp"
#include "deslaurierDubuc.hpp"
#include "waveletMapper.hpp"
    
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
    gp << "set multiplot\n";

    Interval<float> unitInterval(0.0f,1.0f);
    Interval<float> halfInterval(0.25f,0.75f);


    FunctionSample<1000u, float> sample(unitInterval, F);
    //sample.plot(gp);
    
    BinaryTreeNode<float> *tree = new BinaryTreeNode<float>(1,1,unitInterval,nullptr);
    for (unsigned int i = 0; i < 30; i++) {
        tree->insert(sample[i],1);
    } 
    tree->plot(gp);

    tree->computePlacementCondition(1u,0.5);
    tree->plotValid(gp);


    DeslaurierDubuc<float> db(1);
    WaveletMapper<float> simpleDDTreeMapper = [&db](unsigned int j, int k)->Wavelet<float>& { return db; };
    //gp << "clear\n";
    //db.plot(gp,100,1,1);
    
    unsigned int n = tree->countValidNodes();
    std::cout << n << std::endl;
    
    using Eigen::VectorXf;
    using Eigen::MatrixXf;
    MatrixXf A = MatrixXf::Random(n,n);
    VectorXf b = VectorXf::Random(n);
    tree->fillSystem(A,b,simpleDDTreeMapper);


    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    std::cout << A.format(CleanFmt) << std::endl;

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
