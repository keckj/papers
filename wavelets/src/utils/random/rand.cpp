
#include "rand.hpp"
#include <iostream>

namespace Random {
            
                void init() {
                    std::cout << "[Random] Seed init !" << std::endl;
                    srand(time(NULL));
                }

                float randf() {
                        return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                }

                float randf(float LO, float HI) {
                        return LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
                }
                
               int randi(int LO, int HI) {
                        return LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
                }
};
