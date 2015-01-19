#ifndef VEC2_H
#define VEC2_H

#include "headers.hpp"
#include "vec.hpp"

/*
* 2D vector structure of arithmetic type T
* See vec.hpp for more informations and restrictions on T
*/

template <typename T> 
struct Vec2 final : public Vec<2u,T> {
    T& x = this->data[0];
    T& y = this->data[1];

    Vec2();
    Vec2(const Vec2<T> &v);
    Vec2(const Vec<2u,T> &v);
    explicit Vec2(const T data[]);
    explicit Vec2(T x, T y);
    ~Vec2();

    Vec2<T>& operator= (const Vec2<T> &other);

    void setValue(T x, T y);
};
    

template <typename T>
Vec2<T>::Vec2() : Vec<2u,T>() {}

template <typename T>
Vec2<T>::Vec2(T x, T y) : Vec<2u,T>() {
    this->x = x;
    this->y = y;
}
    
template <typename T>
Vec2<T>::Vec2(const Vec2<T> &v) : Vec<2u,T>(v) {
}

template <typename T>
Vec2<T>::Vec2(const Vec<2u,T> &v) : Vec<2u,T>(v) {
}

template <typename T>
Vec2<T>::Vec2(const T data[]) : Vec<2u,T>(data) {}

template <typename T>
Vec2<T>::~Vec2() {}

template <typename T>
Vec2<T> & Vec2<T>::operator= (const Vec2<T> &other) {
    for (unsigned int i = 0; i < 2u; i++) {
        this->data[i] = other.data[i];
    }
    return *this;
}

template <typename T>
void Vec2<T>::setValue(T x, T y) {
    this->x = x;
    this->y = y;
}

#endif /* end of include guard: VEC2_H */
