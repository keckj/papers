
#include "headers.hpp"

#include <iostream>
#include <Eigen/Dense>

#include "rand.hpp"
#include "interval.hpp"
#include "integerGrid.hpp"
#include "binaryTreeNode.hpp"
#include "sample.hpp"
#include "functionSample.hpp"
#include "randomSample.hpp"
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
    return cos(75*t)*exp(t)*cos(10*t);
}

int main(int argc, char **argv) {

    init();
    header();
  
    
    // configuration
    constexpr unsigned int nData = 200u;
    constexpr unsigned int nDataSample = 10000u;
    const Interval<float> interval(0.0f,1.0f);

    //generate and plot samples
    FunctionSample<nDataSample, float> samplePlot(interval, F);
    RandomSample<nData, float> sample(interval, F);

    const PlotBox<float>& box = PlotBox<float>(interval.inf,interval.sup,samplePlot.min(),samplePlot.max()).dilate(1.0);
    //build tree with samples (ALLOCATION)
    TreeNode<float> *tree = new IntegerGrid<float>(interval);

    for (unsigned int i = 0; i < nData; i++) {
        tree->insert(sample[i]);
    } 
    
    //remove bad nodes (check GRP criterion)
    tree->computePlacementCondition(1,0.5);
    
    //map wavelets to tree nodes
    DeslaurierDubuc<float,10> db;
    WaveletMapper<float> simpleDDTreeMapper = [&db](unsigned int j, int k)->Wavelet<float>& { return db; };
    unsigned int n = tree->countValidNodes();
    
    std::cout << n << " samples out of " << nData << " met the GRP criterion !" << std::endl;

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

    FunctionSample<nDataSample, float> reconstructedSample(interval, *waveletTree);
    float L2_distance = (samplePlot - reconstructedSample).norm();
    std::cout << "l2 norm ||F' - F|| = " << L2_distance << " !" << std::endl;

    //plot everything
    Gnuplot gp("tee plot.gp | gnuplot -persist");
    gp << "set term wxt dashed\n";
    gp << "set multiplot\n";
    samplePlot.plotLine(gp, box);
    sample.plotPoints(gp, box);
    tree->plotValidPoints(gp, box);
    waveletTree->plot(gp, box, nDataSample);

    gp << "unset multiplot\n";
    gp << "set term wxt 1 dashed\n";
    gp << "set multiplot\n";
    gp << "unset label\n";
    tree->plot(gp, box, true, true);
    tree->plotValid(gp, box);

    //gp << "unset multiplot\n";
    //gp << "set term wxt 2\n";
    //gp << "set multiplot\n";
    //gp << "unset label\n";
    //waveletTree->plot(gp, box, nDataSample);

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
    
   
