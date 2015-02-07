
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
#include "plotUtils.hpp"
#include "deslaurierDubuc.hpp"
#include "waveletMapper.hpp"
#include "waveletTree.hpp"
#include "plotBox.hpp"
#include <time.h>

void init();
void header();
void footer();
template <unsigned int N, unsigned int P, typename T> 
    T bench();

//interpolated function
float F(float t) {
    return cos(75*t)*exp(t)*cos(10*t);
}
__attribute__((optimize("unroll-loops")))
int main(int argc, char **argv) {

    init();
    header();
    
    DeslaurierDubuc<float, 1u> db1;
    bench<10, 1, float>();
    bench<50, 1, float>();
    bench<100, 1, float>();
    bench<200, 1, float>();
    bench<500, 1, float>();
   
    std::cout << std::endl;
    DeslaurierDubuc<float, 2u> db2;
    bench<10, 2, float>();
    bench<50, 2, float>();
    bench<100, 2, float>();
    bench<200, 2, float>();
    bench<500, 2, float>();
    
    std::cout << std::endl;
    DeslaurierDubuc<float, 3u> db3;
    bench<10, 3, float>();
    bench<50, 3, float>();
    bench<100, 3, float>();
    bench<200, 3, float>();
    bench<500, 3, float>();
    
    std::cout << std::endl;
    DeslaurierDubuc<float, 5u> db5;
    bench<10, 5, float>();
    bench<50, 5, float>();
    bench<100, 5, float>();
    bench<200, 5, float>();
    bench<500, 5, float>();
    
    std::cout << std::endl;
    DeslaurierDubuc<float, 10u> db10;
    bench<10, 10, float>();
    bench<50, 10, float>();
    bench<100, 10, float>();
    bench<200, 10, float>();
    bench<500, 10, float>();
    
    std::cout << std::endl;
    DeslaurierDubuc<float, 50u> db50;
    bench<10, 50, float>();
    bench<50, 50, float>();
    bench<100, 50, float>();
    bench<200, 50, float>();
    bench<500, 50, float>();

    footer();
    return EXIT_SUCCESS;
}


    void header() {
        std::cout << "Bench:" << std::endl << std::endl;
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

    template <unsigned int N, unsigned int P, typename T> 
    T bench() {

        std::cout << "Benching " << N << " random samples with order " << P << "... ";

        double meanError=0, meanInfError=0;
        const unsigned int nReconstructions = 1000u;
            
        // configuration
        constexpr unsigned int nDataSample = 10000u;
        constexpr unsigned int nData = N;
        const Interval<T> interval(0.0f,1.0f);
            
        FunctionSample<nDataSample, T> samplePlot(interval, F);

        clock_t tStart = clock();

        // compute mean error
        for (unsigned int i = 0; i < nReconstructions; i++) {
            //generate random samples
            RandomSample<nData, T> sample(interval, F);

            //build tree with samples (ALLOCATION)
            TreeNode<T> *tree = new IntegerGrid<T>(interval);

            for (unsigned int i = 0; i < nData; i++) {
                tree->insert(sample[i]);
            } 

            //remove bad nodes (check GRP criterion)
            tree->computePlacementCondition(1,0.5);

            //map wavelets to tree nodes
            DeslaurierDubuc<T,P> db;
            WaveletMapper<T> simpleDDTreeMapper = [&db](unsigned int j, int k)->Wavelet<T>& { return db; };
            unsigned int n = tree->countValidNodes();

            //compute wavelet coefficients
            using Eigen::VectorXf;
            using Eigen::MatrixXf;
            MatrixXf A = MatrixXf::Random(n,n);
            VectorXf b = VectorXf::Random(n);

            tree->fillSystem(A,b,simpleDDTreeMapper);
            VectorXf coefficients = A.partialPivLu().solve(b);

            //build wavelet tree with the coefficients
            WaveletTree<T> *waveletTree = WaveletTree<T>::makeTree(tree, simpleDDTreeMapper, coefficients);

            FunctionSample<nDataSample, T> reconstructedSample(interval, *waveletTree);
            double error = (samplePlot - reconstructedSample).norm();
            double infError = (samplePlot - reconstructedSample).normInf();

            if(std::isnan(infError) || !std::isfinite(infError) || infError > 10000
                    || std::isnan(error) || !std::isfinite(error))
                i--;
            else {
                if(i % (nReconstructions/10) == 0)
                    std::cout << PlotUtils::prettyValue<T,3>(T(i)/nReconstructions*100) << "%  " << std::flush;
                meanError += error;
                meanInfError += infError;
            }
        }

        clock_t tEnd = clock();
        double meanTime = static_cast<double>(tEnd - tStart)*1000/(nReconstructions*CLOCKS_PER_SEC);

        std::cout << "100%";
        std::cout << "  => " << meanError/nReconstructions 
            << " : " << meanError/(nReconstructions*nDataSample)
            << " : " << meanInfError/nReconstructions
            << " : " << PlotUtils::prettyValue<double,0>(meanTime) << "ms" << std::endl;;
        return meanError/nReconstructions;
    }


