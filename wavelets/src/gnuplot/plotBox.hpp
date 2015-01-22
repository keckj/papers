
#ifndef PLOTBOX_H
#define PLOTBOX_H

#include "point.hpp"
#include "gnuplot.hpp"

template <typename T>
struct PlotBox {

    Point<T> _xmin;
    Point<T> _xmax;

    PlotBox(Point<T> xmin, Point<T> xmax) :
        _xmin(xmin), _xmax(xmax) {
            assert(_xmin.x < _xmax.x && _xmin.y < _xmax.y);
        }

    PlotBox(T xmin_, T xmax_, T ymin_, T ymax_) :
        _xmin(xmin_, ymin_), _xmax(xmax_, ymax_) {
            assert(xmin_ < xmax_ && ymin_ < ymax_);
    }

    T xmin() const { return _xmin.x; }
    T xmax() const { return _xmax.x; }
    T ymin() const { return _xmin.y; }
    T ymax() const { return _xmax.y; } 

    T width() const { return _xmax.x - _xmin.x; }
    T height() const { return _xmax.y - _xmin.y; }
};

template <typename T>
Gnuplot& operator<< (Gnuplot& gp, const PlotBox<T> &box) {
    gp << "set xrange [" << box.xmin() << ":" << box.xmax() << "]\n";
    gp << "set yrange [" << box.ymin() << ":" << box.ymax() << "]\n";
    return gp;
}

#endif /* end of include guard: PLOTBOX_H */
