
#ifndef SAMPLE_H
#define SAMPLE_H

#include "interval.hpp"
#include "point.hpp"
#include "gnuplot.hpp"
#include <vector>

template <unsigned int N, typename T>
struct Sample {
   
    explicit Sample(const Interval<T> &interval, T x[] = nullptr, T y[] = nullptr);
    explicit Sample(const Interval<T> &interval, std::function<Point<T>(unsigned int, Interval<T>)> lambda);
    Sample(const Sample<N,T> &other);
    virtual ~Sample();

    Interval<T> interval;

    Point<T> data[N];

    void operator+=(const Sample<N,T> &b);
    void operator+=(T b);

    void operator*=(const Sample<N,T> &b);
    void operator*=(T b);

    Point<T> operator[](unsigned int k);

    void plot(Gnuplot &gp) const;
};
    
template <unsigned int N, typename T>
Sample<N,T>::Sample(const Interval<T> &interval, std::function<Point<T>(unsigned int, Interval<T>)> lambda) :
    interval(interval)
{
        for (unsigned int i = 0; i < N; i++) {
            Point<T> sample = lambda(i, interval);
            this->data[i] = sample;
        }
}

template <unsigned int N, typename T>
Sample<N,T>::Sample(const Interval<T> &interval, T _x[], T _y[]) :
    interval(interval) {
        for (unsigned int i = 0; i < N; i++) {
            this->data[i].x = (_x == nullptr ? T(0) : _x[i]);
            this->data[i].y = (_y == nullptr ? T(0) : _y[i]);
        }
}

template <unsigned int N, typename T>
Sample<N,T>::Sample(const Sample<N,T> &other) :
    interval(other.interval) {
}

template <unsigned int N, typename T>
Sample<N,T>::~Sample() {
}

template <unsigned int N, typename T>
void Sample<N,T>::operator+=(const Sample<N,T> &b) {
    for(unsigned int i = 0; i < N; i++) {
        data[i].y += b.data[i].y;
    }
}

template <unsigned int N, typename T>
void Sample<N,T>::operator+=(T b) {
    for(unsigned int i = 0; i < N; i++) {
        data[i].y += b;
    }
}

template <unsigned int N, typename T>
void Sample<N,T>::operator*=(const Sample<N,T> &b) {
    for(unsigned int i = 0; i < N; i++) {
        data[i].y *= b.data[i].y;
    }
}

template <unsigned int N, typename T>
void Sample<N,T>::operator*=(T b) {
    for(unsigned int i = 0; i < N; i++) {
        data[i].y *= b;
    }
}

template <unsigned int N, typename T>
Sample<N,T> operator+(const Sample<N,T>&a, const Sample<N,T> &b) {
    Sample<N,T> sample(a);
    for(unsigned int i = 0; i < N; i++) {
        sample.data[i].x = a.data[i].x;
        sample.data[i].y = a.data[i].y + b.data[i].y;
    }
    return sample;
}

template <unsigned int N, typename T>
Sample<N,T> operator*(const Sample<N,T>&a, const Sample<N,T> &b) {
    Sample<N,T> sample(a);
    for(unsigned int i = 0; i < N; i++) {
        sample.data[i].x = a.data[i].x;
        sample.data[i].y = a.data[i].y * b.data[i].y;
    }
    return sample;
}

template <unsigned int N, typename T>
Sample<N,T> operator*(T a, const Sample<N,T> &b) {
    Sample<N,T> sample(a);
    for(unsigned int i = 0; i < N; i++) {
        sample.data[i].x = a.data[i].x;
        sample.data[i].y = a * b.data[i].y;
    }
    return sample;
}

template <unsigned int N, typename T>
Sample<N,T> operator*(const Sample<N,T>&a, T b) {
    Sample<N,T> sample(a);
    for(unsigned int i = 0; i < N; i++) {
        sample.data[i].x = a.data[i].x;
        sample.data[i].y = a.data[i].y * b;
    }
    return sample;
}
    
template <unsigned int N, typename T>
Point<T> Sample<N,T>::operator[](unsigned int k) {
    assert(k < N);
    return Point<T>(this->data[k]);
}

template <unsigned int N, typename T>
void Sample<N,T>::plot(Gnuplot &gp) const {

    std::vector<std::tuple<float, float>> pts;

    Point<T> Xmin, Xmax;
    for (unsigned int i = 0; i < N; i++) {
        Point<T> X = this->data[i];
        Xmin.x = (X.x < Xmin.x ? X.x : Xmin.x);
        Xmin.y = (X.y < Xmin.y ? X.y : Xmin.y);
        Xmax.x = (X.x > Xmax.x ? X.x : Xmax.x);
        Xmax.y = (X.y > Xmax.y ? X.y : Xmax.y);
        pts.push_back(std::make_tuple(X.x,X.y));
    }

    gp << "set xrange [" << interval.inf << ":" << interval.sup << "]\n";
    gp << "set yrange [" << Xmin.y << ":" << Xmax.y << "]\n";
    gp << "plot '-' title 'samples'\n";
    gp.send1d(pts);
}


#endif /* end of include guard: SAMPLE_H */
