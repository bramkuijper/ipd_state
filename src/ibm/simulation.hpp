#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include <vector>
#include <functional>
#include <random>
#include <iostream>
#include <fstream>
#include "individual.hpp"
#include "parameters.hpp"

class Simulation
{
    private:
        // a data file containing the results of the simulation
        std::ofstream data_file;

        // keep track of the time step of the simulation
        long unsigned time_step = 0;
        
        // random device which is used to generate
        // proper random seeds
        std::random_device rd;

        // store the random seed
        // we need to store this so that we can output the 
        // random seed, so that we could 'replay' the exact
        // same sequence of random numbers for debugging purposes etc
        unsigned int seed;
        
        // random number generator
        std::mt19937 rng_r;
        
        // parameter object
        Parameters params;

        // uniform distribution to compare against probabilities
        std::uniform_real_distribution<double> uniform;
        // uniform dist between + / -
        std::uniform_real_distribution<double> uniform_pm;

        std::vector<Individual> singles;
        std::vector<Individual> offspring;
        std::vector<Individual> paired;

        int number_mortalities = 0;

        void pair_up();
        void interact();
        void reproduce();
        double payoff_pd(double const x, double const xprime);
        double payoff_snowdrift(double const x, double const xprime);
        double payoff(double const x, double const xprime);

        void dismiss_partner();

        int calculate_mortalities();

        void mortality();

        void write_data();
        void write_data_headers();
        void write_parameters();

        void check_state();

    public:
        Simulation(Parameters const &params);

        void run();

};

#endif
