
#ifndef DELAURIERDUBUC_H
#define DELAURIERDUBUC_H

#include <cassert>
#include "wavelet.hpp"
#include "gnuplot.hpp"

template <typename T>
class DeslaurierDubuc : public Wavelet<T> {

    public:
        explicit DeslaurierDubuc(unsigned int order);
        DeslaurierDubuc(const DeslaurierDubuc<T> &other);
        ~DeslaurierDubuc();

        T operator()(T x) const override;
        T operator()(unsigned int j, int k, T x) const override;

        unsigned int order() const;

        static constexpr Interval<T> computeSupport(unsigned int order) {
            return Interval<T>(-T(2)*order + 1/T(2) ,T(2)*order - 1/T(2));
        }

    protected:
        unsigned int _order;
};
        
template <typename T>
DeslaurierDubuc<T>::DeslaurierDubuc(unsigned int order) : 
    Wavelet<T>(DeslaurierDubuc<T>::computeSupport(order)), _order(order) {
        assert(order >= 1);
}

template <typename T>
DeslaurierDubuc<T>::DeslaurierDubuc(const DeslaurierDubuc<T> &other) : Wavelet<T>(other), _order(other.order()) {
}

template <typename T>
DeslaurierDubuc<T>::~DeslaurierDubuc() {
}

template <typename T>
unsigned int DeslaurierDubuc<T>::order() const {
    return _order;
}

template <typename T>
T DeslaurierDubuc<T>::operator()(T x) const {
    if(this->support().contains(x)) {
        return T(1) - sqrt(x*x);
    }
    else {
        return T(0);
    }
}

template <typename T>
T DeslaurierDubuc<T>::operator()(unsigned int j, int k, T x) const {
    return this->operator()(T(std::pow(2,j))*x - T(k));
}

#endif /* end of include guard: DESLAURIERDUBUC_H */
