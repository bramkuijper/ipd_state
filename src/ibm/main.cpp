#include <string>
#include "simulation.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{
    Parameters params;
    params.max_time = std::stoi(argv[1]);

    Simulation sim(params);

    sim.run();
}
