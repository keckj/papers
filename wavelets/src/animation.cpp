
#include "headers.hpp"

#include <iostream>
#include <Eigen/Dense>
#include <boost/filesystem.hpp>
#include <Magick++.h> 
#include <list>
#include <chrono>
#include <thread>

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

    Magick::InitializeMagick(*argv);

    init();
    header();
    
    // configuration
    constexpr unsigned int maxData = 200u;              // 100 samples
    constexpr unsigned int addedPointsPerFrame = 1u;    // 1 additional sample per frame
    constexpr unsigned int maxDeltaT = 25u;             // 50 milisec. animation dt for first frames
    constexpr unsigned int minDeltaT = 1u;              // 10 milisec. animation dt for last frame
    constexpr double power = 20.0;                      // dt decreases with cos(x*pi/2)^power
    constexpr unsigned int order = 10u;                 // order of the DD wavelet

    constexpr unsigned int nDataSample = 10000u;
    const Interval<float> interval(0.0f,1.0f);
   
    //others
    constexpr unsigned int nFrames = ceil(maxData/addedPointsPerFrame);

    //generate and plot samples
    FunctionSample<nDataSample, float> samplePlot(interval, F);
    RandomSample<maxData+addedPointsPerFrame, float> sample(interval, F);

    const PlotBox<float>& box = PlotBox<float>(interval.inf,interval.sup,samplePlot.min(),samplePlot.max()).dilate(1.0);
    DeslaurierDubuc<float,order> db;
        
    //initialize gnuplot for png output
    Gnuplot gp("tee animate.gp | gnuplot -persist 2>/dev/null");
    gp << "set terminal pngcairo size 800,600 enhanced font 'Helvetica,20' dashed\n";
    gp << "unset label\n";

    //create folder
	boost::filesystem::path dir("img");
	assert(boost::filesystem::is_directory(dir) || boost::filesystem::create_directories(dir));

    //image container
    std::list<Magick::Image> images;
    Magick::Image img;
        
    std::cout << "Generating frames :" << std::endl;
    // GNUPLOT CORRUPT LAST GENERATED IMAGE -- so just create one more
    for (unsigned int i = 1; i <= nFrames+1; i++) {
        if(i == nFrames + 1)
            std::cout << "\tGenerating dummy frame " << i << "..." << std::endl;
        else
            std::cout << "\tGenerating frame " << i << "..." << std::endl;
        
        unsigned int nData = i*addedPointsPerFrame;
        if(i == nFrames)
            nData = maxData;

        //build tree with samples (ALLOCATION)
        TreeNode<float> *tree = new IntegerGrid<float>(interval);
        for (unsigned int k = 0; k < nData; k++) {
            tree->insert(sample[k]);
        } 

        //remove bad nodes (check GRP criterion)
        tree->computePlacementCondition(1,0.5);

        //map wavelets to tree nodes
        WaveletMapper<float> simpleDDTreeMapper = [&db](unsigned int j, int k)->Wavelet<float>& { return db; };
        unsigned int n = tree->countValidNodes();

        //compute wavelet coefficients
        using Eigen::VectorXf;
        using Eigen::MatrixXf;
        MatrixXf A = MatrixXf::Random(n,n);
        VectorXf b = VectorXf::Random(n);

        tree->fillSystem(A,b,simpleDDTreeMapper);
        VectorXf coefficients = A.partialPivLu().solve(b);

        //build wavelet tree with the coefficients
        WaveletTree<float> *waveletTree = WaveletTree<float>::makeTree(tree, simpleDDTreeMapper, coefficients);

        FunctionSample<nDataSample, float> reconstructedSample(interval, *waveletTree);
   
        std::string imgPath = "img/interpolation_" + std::to_string(i) + ".png";
        gp << "unset multiplot\n";
        gp << "set output '" << imgPath << "'\n";
        gp << "set multiplot\n";
        samplePlot.plotLine(gp, box);
        sample.plotPoints(gp, box, nData);
        tree->plotValidPoints(gp, box);
        waveletTree->plot(gp, box, nDataSample);
    }

    std::cout << std::endl;
    std::cout << "Generating animation..." << std::endl;
    
    for (unsigned int i = 1; i <= nFrames; i++) {
        std::string imgPath = "img/interpolation_" + std::to_string(i) + ".png";
        img.read(imgPath);
        double x = static_cast<double>(i)/(nFrames); 
        unsigned int dT = static_cast<unsigned int>(maxDeltaT * std::pow(cos(x*3.14159265/2),power));
        dT = std::max(dT, minDeltaT);
        img.animationDelay(dT);
        images.push_back(img);
    }

    Magick::writeImages(images.begin(), images.end(), "img/interpolation.gif");

    std::cout << "Animation generated at img/interpolation.gif !" << std::endl;

    footer();

    return EXIT_SUCCESS;
}
    
void header() {
    std::cout << std::endl;
    std::cout << "Generating animation !" << std::endl;
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
    
   
