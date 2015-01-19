
#ifndef FUNCTIONSAMPLE_H
#define FUNCTIONSAMPLE_H

#include "rand.hpp"

template <unsigned int N, typename T>
struct FunctionSample : public Sample<N,T> {
    explicit FunctionSample(const Interval<T> &interval, std::function<T(T)> lambda);
    FunctionSample(const Sample<N,T> &other);
    virtual ~FunctionSample();
};
    
template <unsigned int N, typename T>
FunctionSample<N,T>::FunctionSample(const Interval<T> &interval, std::function<T(T)> lambda) :
    Sample<N,T>(interval) {
        for (unsigned int i = 0; i < N; i++) {
            T x = interval.inf + Random::randf()*(interval.sup - interval.inf);
            T y = lambda(x);
            this->data[i] = Point<T>(x,y);
        }
}

template <unsigned int N, typename T>
FunctionSample<N,T>::FunctionSample(const Sample<N,T> &other) :
    Sample<N,T>(other) {
}

template <unsigned int N, typename T>
FunctionSample<N,T>::~FunctionSample() {
}

#endif
