
#include "deslaurierDubucUtils.hpp"

#include <boost/math/special_functions/binomial.hpp>

namespace DeslaurierDubucUtils {
        unsigned long binom(unsigned long n, unsigned long k) {
            return static_cast<unsigned long>(boost::math::binomial_coefficient<double>(n,k));
        }
}
