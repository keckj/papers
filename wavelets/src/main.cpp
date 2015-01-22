
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
#include "waveletTree.hpp"
#include "plotBox.hpp"

void init();
void header();
void footer();

//interpolated function
float F(float t) {
    return exp(t);
}

int main(int argc, char **argv) {

    init();
    header();
  
    
    // configuration
    constexpr unsigned int nData = 100u;
    constexpr unsigned int nDataSample = 1000u;
    Interval<float> interval(0.0f,1.0f);

    
    //generate and plot samples
    FunctionSample<nDataSample, float> samplePlot(interval, F);
    FunctionSample<nData, float> sample(interval, F);

    const PlotBox<float>& box = PlotBox<float>(interval.inf,interval.sup,samplePlot.min(),samplePlot.max());
    Gnuplot gp("tee plot.gp | gnuplot -persist");
    gp << "set multiplot\n";

    samplePlot.plotLine(gp, box);
    sample.plotPoints(gp, box);

    
    //build tree with samples (ALLOCATION)
    TreeNode<float> *tree = new BinaryTreeNode<float>(1,1,interval,nullptr);
    for (unsigned int i = 0; i < nData; i++) {
        tree->insert(sample[i],1);
    } 
    tree->plot(gp, box, true, true);

    return EXIT_SUCCESS;

    //remove bad nodes (check GRP criterion)
    tree->computePlacementCondition(1u,0.5);

    tree->plotValidPoints(gp, box);
    tree->plotValid(gp, box);

    
    //map wavelets to tree nodes
    DeslaurierDubuc<float> db(1);
    WaveletMapper<float> simpleDDTreeMapper = [&db](unsigned int j, int k)->Wavelet<float>& { return db; };
    unsigned int n = tree->countValidNodes();
    
    std::cout << n << " samples out of " << nData << " met the GRP criterion !";
   

    //compute wavelet coefficients
    using Eigen::VectorXf;
    using Eigen::MatrixXf;
    MatrixXf A = MatrixXf::Random(n,n);
    VectorXf b = VectorXf::Random(n);

    tree->fillSystem(A,b,simpleDDTreeMapper);
    VectorXf coefficients = A.partialPivLu().solve(b);
    
    float epsilon = (A*coefficients -b).maxCoeff();
    std::cout << "Maximum sample error : " << epsilon << std::endl;


    //build wavelet tree with the coefficients
    WaveletTree<float> *waveletTree = WaveletTree<float>::makeTree(tree, simpleDDTreeMapper, coefficients);

    gp << "unset multiplot\n";
    gp << "set term wxt 1\n";
    gp << "unset label\n";
    waveletTree->plot(gp, box, nDataSample);


    footer();

    return EXIT_SUCCESS;
}
    
//Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
//std::cout << x.format(CleanFmt) << std::endl;
    

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

void init() {
    Random::init();

    Globals::init();
    Globals::check();
}
    
   