#include <string>
#include "simulation.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{
    Parameters params;
    params.max_time = std::stoi(argv[1]);
    params.init_x = std::stod(argv[2]);
    params.init_y = std::stod(argv[3]);
    params.mu_y = std::stod(argv[4]);
    params.mu_xp = std::stod(argv[5]);
    params.mu_yp = std::stod(argv[6]);
    params.cooperation_type = static_cast<game_type>(std::stoi(argv[7]));
    params.startup_cost = std::stod(argv[8]);
    params.mortality_prob = std::stod(argv[9]);
    params.resource_variation = std::stod(argv[10]);
    params.file_base_name = argv[11];

    Simulation sim(params);

    sim.run();
}
