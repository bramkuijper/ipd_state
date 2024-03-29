#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

enum game_type
{
    prisoners_dilemma=0,
    snowdrift_game=1
};

class Parameters
{
    public:
        // probability that the 
        // interaction will be continued

        double mu_x = 0.02;
        double mu_y = 0.02;
        double mu_xp = 0.01;
        double mu_yp = 0.01;

        game_type cooperation_type = prisoners_dilemma;

        double sdmu = 0.01;

        unsigned int N = 1000;
        long int max_time = 10000;

        double init_x = 0.5;
        double init_y = 0.5;

        double startup_cost = 0.01;

        double mortality_prob = 0.1;

        double dismiss_error = 0.01;

        double resource_variation = 1;

        double cost_of_reproduction = 2;

        int data_interval = 1;
        int data_interval_whole_pop = 100;

        std::string file_base_name = "data_ipd";

        double max_resources = 100.0;

        bool whole_population{false};
};

#endif
