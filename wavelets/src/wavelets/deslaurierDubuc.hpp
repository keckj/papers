
#ifndef DELAURIERDUBUC_H
#define DELAURIERDUBUC_H

#include <cassert>
#include "wavelet.hpp"
#include "gnuplot.hpp"
//#include <unsupported/Eigen/FFT>
#include "deslaurierDubucUtils.hpp"

template <typename T, unsigned int N>
class DeslaurierDubuc : public Wavelet<T> {

    public:
        DeslaurierDubuc();
        DeslaurierDubuc(const DeslaurierDubuc<T,N> &other);
        ~DeslaurierDubuc();

        T operator()(T x) const override;
        T operator()(unsigned int j, int k, T x) const override;

        static constexpr unsigned int sampleCount();

        static std::vector<T> getSamples();

    protected:
        static constexpr unsigned int _unitSamplesCount = 100u;
        static std::vector<T> _samples;
        static bool _init;
        
        static constexpr Interval<T> computeSupport();

        static void init();
        static unsigned int toSampleId(T x);
};

        
template <typename T, unsigned int N>
bool DeslaurierDubuc<T,N>::_init = false;

template <typename T, unsigned int N>
std::vector<T> DeslaurierDubuc<T,N>::_samples;

template <typename T, unsigned int N>
void DeslaurierDubuc<T,N>::init() {
    if(_init)
        return;
}
        
template <typename T, unsigned int N>
unsigned int DeslaurierDubuc<T,N>::toSampleId(T x) {
    Interval<T> sup = computeSupport();
    assert(sup.contains(x));
    return static_cast<unsigned int>((x-sup.inf)/sup.length() * (sampleCount() - 1u));
}
        
template <typename T, unsigned int N>
DeslaurierDubuc<T,N>::DeslaurierDubuc() : 
    Wavelet<T>(DeslaurierDubuc<T,N>::computeSupport()) {
        init();
}

template <typename T, unsigned int N>
DeslaurierDubuc<T,N>::DeslaurierDubuc(const DeslaurierDubuc<T,N> &other) : Wavelet<T>(other) {
}

template <typename T, unsigned int N>
DeslaurierDubuc<T,N>::~DeslaurierDubuc() {
}

template <typename T, unsigned int N>
T DeslaurierDubuc<T,N>::operator()(T x) const {
    if(this->support(0,0).contains(x)) {
        return T(1) - sqrt(x*x);
    }
    else {
        return T(0);
    }
}

template <typename T, unsigned int N>
T DeslaurierDubuc<T,N>::operator()(unsigned int j, int k, T x) const {
    if(N == 0)
        return this->operator()(T(std::pow(2,j))*x - T(k));
    else 
        return 0;
}
        
template <typename T, unsigned int N>
constexpr Interval<T> DeslaurierDubuc<T,N>::computeSupport() {
    switch(N) {
        case(0u): 
            return Interval<T>(-1,1);
        default:
            return Interval<T>(-T(2)*N + 1/T(2) ,T(2)*N - 1/T(2));
    }
}
        
template <typename T, unsigned int N>
std::vector<T> DeslaurierDubuc<T,N>::getSamples() {
    return _samples;
}

template <typename T, unsigned int N>
constexpr unsigned int DeslaurierDubuc<T,N>::sampleCount() {
    return static_cast<unsigned int>(computeSupport().length()*_unitSamplesCount);
}

#endif /* end of include guard: DESLAURIERDUBUC_H */
