
#ifndef DELAURIERDUBUC_H
#define DELAURIERDUBUC_H

#include <cassert>
#include "wavelet.hpp"
#include "gnuplot.hpp"
//#include <unsupported/Eigen/FFT>
#include "deslaurierDubucUtils.hpp"

template <typename T, unsigned int P>
class DeslaurierDubuc : public Wavelet<T> {

    public:
        DeslaurierDubuc();
        DeslaurierDubuc(const DeslaurierDubuc<T,P> &other);
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

        
template <typename T, unsigned int P>
bool DeslaurierDubuc<T,P>::_init = false;

template <typename T, unsigned int P>
std::vector<T> DeslaurierDubuc<T,P>::_samples;

template <typename T, unsigned int P>
T DeslaurierDubuc<T,P>::_resolution = T(0);

template <typename T, unsigned int P>
unsigned int DeslaurierDubuc<T,P>::_levels = 0u;

template <typename T, unsigned int P>
void DeslaurierDubuc<T,P>::init() {
    if(_init)
        return;

    _levels = 14u;
    _resolution = T(1)/std::pow<T>(2, _levels);
    std::cout << "Generaring Deslaurier-Dubuc wavelet of order " << P << " with " << _levels 
        << " convolutions (" << (1u<<_levels) << " samples per unit, resolution dx=" << _resolution << ")..." << std::endl;
    _samples = DeslaurierDubucUtils::generateScalingFunction<T>(P, _levels);
    _init = true;
}
        
template <typename T, unsigned int P>
unsigned int DeslaurierDubuc<T,P>::toSampleId(T x) {
    assert(_init);
    Interval<T> sup = DeslaurierDubucUtils::computeSupport<T>(P);
    assert(sup.contains(x));
    return floor((x-sup.inf)/sup.length() * (sampleCount() - 1u));
}
        
template <typename T, unsigned int P>
DeslaurierDubuc<T,P>::DeslaurierDubuc() : 
    Wavelet<T>(DeslaurierDubucUtils::computeSupport<T>(P)) {
        init();
}

template <typename T, unsigned int P>
DeslaurierDubuc<T,P>::DeslaurierDubuc(const DeslaurierDubuc<T,P> &other) : Wavelet<T>(other) {
}

template <typename T, unsigned int P>
DeslaurierDubuc<T,P>::~DeslaurierDubuc() {
}

template <typename T, unsigned int P>
T DeslaurierDubuc<T,P>::operator()(T x) const {

    //std::cout << this->support() << std::endl;
    //std::cout << "x = " << x << std::endl;

    if(! this->support().contains(x)) {
        return T(0);
    }


    unsigned int sampleIdLow, sampleIdHigh;

    sampleIdLow = toSampleId(x);
    if(sampleIdLow == _samples.size() - 1u)
        sampleIdLow--;
    
    sampleIdHigh = sampleIdLow + 1u;

    T inf = this->support().inf;
    T xlow =  inf + sampleIdLow  * _resolution; 
    T xhigh = inf + sampleIdHigh * _resolution; 
    T ylow =  _samples[sampleIdLow];
    T yhigh = _samples[sampleIdHigh];
    
    //std::cout << "xlow = " << xlow << std::endl;
    //std::cout << "xhigh = " << xhigh << std::endl;
    //std::cout << std::endl;

    T a = (yhigh - ylow)/(xhigh - xlow);
    T b = ylow - a * xlow;

    return (a*x + b);
}

template <typename T, unsigned int P>
T DeslaurierDubuc<T,P>::operator()(unsigned int j, int k, T x) const {
        return this->operator()(T(std::pow(2,j))*x - T(k));
}
        
template <typename T, unsigned int P>
std::vector<T> DeslaurierDubuc<T,P>::getSamples() {
    return _samples;
}

template <typename T, unsigned int P>
constexpr unsigned int DeslaurierDubuc<T,P>::sampleCount() {
    return _samples.size();
}

#endif /* end of include guard: DESLAURIERDUBUC_H */
