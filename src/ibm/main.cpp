#include <string>
#include "simulation.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{
    Parameters params;
    params.w = std::stod(argv[1]);

    Simulation sim(params);
}
