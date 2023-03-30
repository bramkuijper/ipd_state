#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"

class Individual
{
    public:
        double x = 0.0; // cooperativeness baseline trait 
        double xp = 0.0; // cooperativeness plasticity slope
        double y = 0.0; // choosiness baseline threshold
        double yp = 0.0; // choosiness plasticity threshold

        double resources = 0.0;
        
        Individual(double const init_x, double const init_y);
        Individual(Individual const &other);

        // birth constructor
        Individual(Individual const &parent
                ,Parameters const &params
                ,std::mt19937 rng_r);

        void operator=(Individual const &other);
};

#endif
