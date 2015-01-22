
#ifndef AFFINETRANSFORMATION_H
#define AFFINETRANSFORMATION_H

#include "plotBox.hpp"
#include "point.hpp"

template <typename T>
class AffineTransformation {

    public:
        AffineTransformation(const PlotBox<T> &b1, const PlotBox<T> &b2) {
           ax = b2.width()/b1.width(); 
           bx = b2.xmin() - ax*b1.xmin();

           ay = b2.height()/b1.height(); 
           by = b2.ymin() - ay*b1.ymin();
        }

        Point<T> operator()(const Point<T> &X) const {
            return Point<T>(ax*X.x+bx, ay*X.y+by);
        }

    protected:

        T ax, ay;
        T bx, by;
};


#endif /* end of include guard: AFFINETRANSFORMATION_H */
