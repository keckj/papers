
#ifndef RAND_H
#define RAND_H

#include <cstdlib>
#include <time.h>

namespace Random {
    void init();
    float randf();
    float randf(float LO, float HI);
    int randi(int LO, int HI);
};

#endif /* end of include guard: RAND_H */
