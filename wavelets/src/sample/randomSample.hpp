
#ifndef RANDOMSAMPLE_H
#define RANDOMSAMPLE_H

#include "rand.hpp"

template <unsigned int N, typename T>
struct RandomSample : public Sample<N,T> {

    explicit RandomSample(const Interval<T> &interval, const std::function<T(T)> &lambda);
    RandomSample(const Sample<N,T> &other);
    virtual ~RandomSample();

    T min() const;
    T max() const;

    protected:
        T _min, _max;

};
    
template <unsigned int N, typename T>
RandomSample<N,T>::RandomSample(const Interval<T> &interval, const std::function<T(T)> &lambda) :
    Sample<N,T>(interval) {
        _min = lambda(interval.inf);
        _max = lambda(interval.inf);
        for (unsigned int i = 0; i < N; i++) {
            T x = interval.inf + Random::randf() * (interval.sup - interval.inf);
            T y = lambda(x);
            this->data[i] = Point<T>(x,y);

            _min = (y < _min ? y : _min);
            _max = (y > _max ? y : _max);
        }
}

template <unsigned int N, typename T>
RandomSample<N,T>::RandomSample(const Sample<N,T> &other) :
    Sample<N,T>(other), _min(min), _max(max) {
}

template <unsigned int N, typename T>
RandomSample<N,T>::~RandomSample() {
}
    
template <unsigned int N, typename T>
T RandomSample<N,T>::min() const {
    return _min;
}

template <unsigned int N, typename T>
T RandomSample<N,T>::max() const {
    return _max;
}

#endif
