
#ifndef DELAURIERDUBUC_H
#define DELAURIERDUBUC_H

#include <cassert>
#include "wavelet.hpp"
#include "gnuplot.hpp"
//#include <unsupported/Eigen/FFT>
#include "deslaurierDubucUtils.hpp"

template <typename T, unsigned int M>
class DeslaurierDubuc : public Wavelet<T> {

    public:
        DeslaurierDubuc();
        DeslaurierDubuc(const DeslaurierDubuc<T,M> &other);
        ~DeslaurierDubuc();

        T operator()(T x) const override;
        T operator()(unsigned int j, int k, T x) const override;

        static constexpr unsigned int sampleCount();
        static std::vector<T> getSamples();

    protected:
        static unsigned int _levels;
        static T _resolution;
        static std::vector<T> _samples;
        static bool _init;
        
        static void init();
        static unsigned int toSampleId(T x);
};

        
template <typename T, unsigned int M>
bool DeslaurierDubuc<T,M>::_init = false;

template <typename T, unsigned int M>
std::vector<T> DeslaurierDubuc<T,M>::_samples;

template <typename T, unsigned int M>
T DeslaurierDubuc<T,M>::_resolution = T(0);

template <typename T, unsigned int M>
unsigned int DeslaurierDubuc<T,M>::_levels = 0u;

template <typename T, unsigned int M>
void DeslaurierDubuc<T,M>::init() {
    if(_init)
        return;

    _levels = 10u;
    _resolution = T(1)/std::pow<T>(2, _levels);
    std::cout << "Generaring Deslaurier-Dubuc wavelet of order " << M << " with " << _levels 
        << " convolutions (" << (1u<<_levels) << " samples per unit, resolution dx=" << _resolution << ")..." << std::endl;
    _samples = DeslaurierDubucUtils::generateScalingFunction<T>(M, _levels);
    _init = true;
}
        
template <typename T, unsigned int M>
unsigned int DeslaurierDubuc<T,M>::toSampleId(T x) {
    assert(_init);
    Interval<T> sup = DeslaurierDubucUtils::computeSupport<T>(M);
    assert(sup.contains(x));
    return floor((x-sup.inf)/sup.length() * (sampleCount() - 1u));
}
        
template <typename T, unsigned int M>
DeslaurierDubuc<T,M>::DeslaurierDubuc() : 
    Wavelet<T>(DeslaurierDubucUtils::computeSupport<T>(M)) {
        init();
}

template <typename T, unsigned int M>
DeslaurierDubuc<T,M>::DeslaurierDubuc(const DeslaurierDubuc<T,M> &other) : Wavelet<T>(other) {
}

template <typename T, unsigned int M>
DeslaurierDubuc<T,M>::~DeslaurierDubuc() {
}

template <typename T, unsigned int M>
T DeslaurierDubuc<T,M>::operator()(T x) const {

    if(! this->support().contains(x)) {
        return T(0);
    }

    unsigned int sampleIdLow, sampleIdHigh;

    sampleIdLow = toSampleId(x);
    if(sampleIdLow == _samples.size())
        sampleIdLow--;
    
    sampleIdHigh = sampleIdLow + 1u;

    T xlow =  sampleIdLow  * _resolution; 
    T xhigh = sampleIdHigh * _resolution; 
    T ylow =  _samples[sampleIdLow];
    T yhigh = _samples[sampleIdHigh];

    T a = (yhigh - ylow)/(xhigh - xlow);
    T b = ylow - a * xlow;

    return (a*x + b);
}

template <typename T, unsigned int M>
T DeslaurierDubuc<T,M>::operator()(unsigned int j, int k, T x) const {
    if(M == 0)
        return this->operator()(T(std::pow(2,j))*x - T(k));
    else 
        return 0;
}
        
template <typename T, unsigned int M>
std::vector<T> DeslaurierDubuc<T,M>::getSamples() {
    return _samples;
}

template <typename T, unsigned int M>
constexpr unsigned int DeslaurierDubuc<T,M>::sampleCount() {
    return _samples.size();
}

#endif /* end of include guard: DESLAURIERDUBUC_H */
