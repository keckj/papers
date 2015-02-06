

#ifndef DESLAURIERDUBUCUTILS_H
#define DESLAURIERDUBUCUTILS_H

#include <vector>
#include "interval.hpp"

namespace DeslaurierDubucUtils {


    unsigned long binom(unsigned long n, unsigned long k);
   
    template <typename T>
    constexpr Interval<T> computeSupport(unsigned int m) {
        switch(m) {
            case(0u): 
                return Interval<T>(-1,1);
            default:
                return Interval<T>(-T(2)*m + 1/T(2) ,T(2)*T(m) - 1/T(2));
        }
    }

    template <typename T>
        std::vector<T> lagrange(unsigned int p, T x) {

            std::vector<T> lagrangeCoefficients;
            T coef;

            for (unsigned int i = 0; i <= 2*p+1; i++) {

                coef = T(1);
                for (unsigned int j = 0; j <= 2*p+1; j++) {
                    if(i != j) {
                        coef *= (x+T(p)-T(j))/(T(i)-T(j));
                    }
                }

                lagrangeCoefficients.push_back(coef);
            }

            return lagrangeCoefficients;
        }


    template <typename T>
        std::vector<T> deslaurierDubuc(unsigned int m, T x) {

            int mm = static_cast<int>(m);

            std::vector<T> lagrangeCoefficients;
            std::vector<T> deslaurierDubucCoefficients;
            T coef;

            lagrangeCoefficients = DeslaurierDubucUtils::lagrange<T>(m/2u, x);

            for (int k = -2*mm+1; k <= 2*mm-1; k++) {
                if(k % 2 == 0) 
                    coef = (k == 0 ? T(1) : T(0));
                else
                    coef = lagrangeCoefficients[m + (k-1)/2];

                deslaurierDubucCoefficients.push_back(coef);
            }

            return deslaurierDubucCoefficients;
        }


    template <typename T>
        std::vector<T> generateScalingFunction(unsigned int m, unsigned int levels) {
            
            const std::vector<T> deslaurierDubucCoefficients = DeslaurierDubucUtils::deslaurierDubuc<T>(m,T(0.5));

            const unsigned int nCoefficients = deslaurierDubucCoefficients.size();
            const int halfCoef = (nCoefficients - 1)/2u;
            const int coefCenter = (nCoefficients - 1u)/2u;

            const Interval<T> support = DeslaurierDubucUtils::computeSupport<T>(m);

            const unsigned int unitLength = (1u << levels);
            const unsigned int length = static_cast<unsigned int>(support.length());
            const unsigned int size =  length * unitLength + 1u;

            //std::cout << support << std::endl;
            //std::cout << "unitlength " << unitLength << std::endl;
            //std::cout << "length " << length << std::endl;
            //std::cout << "size " << size << std::endl;
            
            //init
            std::vector<T> scalingFunction1(size, T(0));
            std::vector<T> scalingFunction2(size);
            scalingFunction1[(size-1u)/2u] = T(1);

            std::vector<T> scalingFunctions[2];
            scalingFunctions[0] = scalingFunction1;
            scalingFunctions[1] = scalingFunction2;

            //levels convolutions 
            unsigned int offset, count;
            int id, iid;
            for(unsigned int level = 1u; level <= levels; level++) {

                std::vector<T> &src = scalingFunctions[(level+1u)%2];
                std::vector<T> &dst = scalingFunctions[level%2];
                offset = (1u << (levels - level));
                count =  (1u << level) * length + 1u;

                for (unsigned int i = 0; i < count; i++) {
                    id = i*offset;
                    dst[id] = T(0);
                    for (int k = -halfCoef; k <= halfCoef; k++) {
                        iid = id + k*offset; 
                        if(iid >= 0 && iid < size)
                            dst[id] += deslaurierDubucCoefficients[coefCenter + k] * src[iid];
                    }
                }
                
                //std::cout << "current level is " << level << "/" << levels << std::endl;
                //std::cout << "offset is " << offset << std::endl;
                //std::cout << "count is " << count << std::endl;
            }

            return scalingFunctions[levels%2];
        }
}


#endif /* end of include guard: DESLAURIERDUBUCUTILS_H */
