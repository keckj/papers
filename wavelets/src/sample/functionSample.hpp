
#ifndef FUNCTIONSAMPLE_H
#define FUNCTIONSAMPLE_H

#include "rand.hpp"

template <unsigned int N, typename T>
struct FunctionSample : public Sample<N,T> {

    explicit FunctionSample(const Interval<T> &interval, const std::function<T(T)> &lambda);
    FunctionSample(const Sample<N,T> &other);
    virtual ~FunctionSample();

    T min() const;
    T max() const;

    protected:
        T _min, _max;

};
    
template <unsigned int N, typename T>
FunctionSample<N,T>::FunctionSample(const Interval<T> &interval, const std::function<T(T)> &lambda) :
    Sample<N,T>(interval) {
        _min = lambda(interval.inf);
        _max = lambda(interval.inf);
        for (unsigned int i = 0; i < N; i++) {
            T x = interval.inf + i*(interval.sup - interval.inf)/T(N-1);
            T y = lambda(x);
            this->data[i] = Point<T>(x,y);

            _min = (y < _min ? y : _min);
            _max = (y > _max ? y : _max);
        }
}

template <unsigned int N, typename T>
FunctionSample<N,T>::FunctionSample(const Sample<N,T> &other) :
    Sample<N,T>(other), _min(min), _max(max) {
}

template <unsigned int N, typename T>
FunctionSample<N,T>::~FunctionSample() {
}
    
template <unsigned int N, typename T>
T FunctionSample<N,T>::min() const {
    return _min;
}

template <unsigned int N, typename T>
T FunctionSample<N,T>::max() const {
    return _max;
}

#endif
