

#ifndef DESLAURIERDUBUCUTILS_H
#define DESLAURIERDUBUCUTILS_H

#include <vector>

namespace DeslaurierDubucUtils {

        unsigned long binom(unsigned long n, unsigned long k);

        template <typename T>
        std::vector<T> lagrange(unsigned int p, T x);
        
        template <typename T>
        std::vector<T> deslaurierDubuc(unsigned int m, T x);
        
        template <typename T>
        std::vector<T> generateScalingFunction(unsigned int m, unsigned int levels);
}

#endif /* end of include guard: DESLAURIERDUBUCUTILS_H */
