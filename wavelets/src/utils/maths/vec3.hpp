
#ifndef VEC3_H
#define VEC3_H

#include "headers.hpp"
#include "vec.hpp"

/*
* 3D vector structure of arithmetic type T
* See vec.hpp for more informations and restrictions on T
*/

template <typename T> 
struct Vec3 : public Vec<3u,T> {
    T& x = this->data[0];
    T& y = this->data[1];
    T& z = this->data[2];

    Vec3();
    Vec3(const Vec3<T> &v);       //mandatory for references buqs
    Vec3(const Vec<3u,T> &v);     //mandatory for operations on Vec<3u,T>
    explicit Vec3(const T data[]);
    explicit Vec3(T x, T y, T z);
    ~Vec3();

    void setValue(T x, T y, T z);
    
    Vec3<T> & operator= (const Vec3<T> &other);

    Vec3<T> & operator^= (const Vec3<T> &a);
    
    Vec3<T> orthogonalVec() const;
};
    

template <typename T>
Vec3<T>::Vec3() : Vec<3u,T>() {
}

template <typename T>
Vec3<T>::Vec3(T x, T y, T z) : Vec<3u,T>() {
    this->x = x;
    this->y = y;
    this->z = z;
}

template <typename T>
Vec3<T>::Vec3(const Vec<3u,T> &v) : Vec<3u,T>(v) {
}
    
template <typename T>
Vec3<T>::Vec3(const Vec3<T> &v) : Vec<3u,T>(v) {
}

template <typename T>
Vec3<T>::Vec3(const T data[]) : Vec<3u,T>(data) {
}

template <typename T>
Vec3<T>::~Vec3() {}

template <typename T>
void Vec3<T>::setValue(T x, T y, T z) {
    this->x = x;
    this->y = y;
    this->z = z;
}
    
template <typename T>
Vec3<T> & Vec3<T>::operator= (const Vec3<T> &other) {
    for (unsigned int i = 0; i < 3u; i++) {
        this->data[i] = other.data[i];
    }
    return *this;
}

template <typename T>
Vec3<T> & Vec3<T>::operator^= (const Vec3<T> &a) {
    Vec3<T> b(*this);
    x = b.y*a.z - b.z*a.y;
    y = b.z*a.x - b.x*a.z;
    z = b.x*a.y - b.y*a.x;
    return *this;
}

template <typename T>
Vec3<T> operator^ (const Vec3<T> &a, const Vec3<T> &b) {
    return Vec3<T>(
            a.y*b.z - a.z*b.y,
            a.z*b.x - a.x*b.z,
            a.x*b.y - a.y*b.x
            );
}

template <typename T>
Vec3<T> Vec3<T>::orthogonalVec () const {
    return Vec3<T>(z,z,-x-y);
}

#endif /* end of include guard: VEC3_H */
