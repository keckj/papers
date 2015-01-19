
#ifndef INTEGERGRID_H
#define INTEGERGRID_H

template <typename T, typename E>
class IntegerGrid {

    public:
        IntegerGrid(const Interval<T> &interval);
        IntegerGrid(const IntegerGrid<T,E> &other);
        IntegerGrid<T,E>& operator= (const IntegerGrid<T,E> &other);
        ~IntegerGrid();
};

template <typename T, typename E>
IntegerGrid<T,E>::IntegerGrid(const Interval<T> &interval) {}

template <typename T, typename E>
IntegerGrid<T,E>::IntegerGrid(const IntegerGrid<T,E> &other) {}

template <typename T, typename E>
IntegerGrid<T,E>& IntegerGrid<T,E>::operator= (const IntegerGrid<T,E> &other) { return *this; }

template <typename T, typename E>
IntegerGrid<T,E>::~IntegerGrid() {}

#endif
