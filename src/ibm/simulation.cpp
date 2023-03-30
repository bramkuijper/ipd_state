#include <iostream>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include "simulation.hpp"
#include "parameters.hpp"


Simulation::Simulation(Parameters const &params) :
    rd{} // initialize random device, see *.h file
    ,seed{rd()} // initialize seed
    ,rng_r{seed} // initialize the random number generator
    ,uniform{0.0,1.0} // initialize the uniform distribution
    ,params{params}
{}
