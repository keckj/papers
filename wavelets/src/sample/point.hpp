
#ifndef POINT_H
#define POINT_H

#include <iostream>

template <typename T>
struct Point {
    T x;
    T y;

    Point() : x(0), y(0) {};
    Point(T x, T y) : x(x), y(y) {};
    Point(const Point<T> &other) : x(other.x), y(other.y) {};
    ~Point() {};
};

template <typename T>
std::ostream & operator<< (std::ostream &os, const Point<T> &point) {
    os << "(" << point.x << "," << point.y << ")";
    return os;
}



#endif /* end of include guard: POINT_H */
