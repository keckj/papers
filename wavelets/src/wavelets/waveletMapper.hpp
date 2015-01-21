
#ifndef WAVELETMAPPER_H
#define WAVELETMAPPER_H

#include <functional>
#include "wavelet.hpp"

template <typename T>
using WaveletMapper = std::function<Wavelet<T>&(unsigned int j, int k)>;

#endif /* end of include guard: WAVELETMAPPER_H */
