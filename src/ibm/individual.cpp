#include <random>
#include "individual.hpp"
#include "parameters.hpp"

Individual::Individual(double const init_x, double const init_y) :
    x{init_x}
    ,y{init_y}
    ,xp{0.0}
    ,yp{0.0}
    ,resources{0.0}
{}

Individual::Individual(Individual const &other) :
    x{other.x}
    ,y{other.y}
    ,xp{other.xp}
    ,yp{other.yp}
    ,resources{other.resources}
{
}

Individual::Individual(Individual const &parent
        ,Parameters const &params
        ,std::mt19937 rng_r) :
    x{parent.x}
    ,y{parent.y}
    ,xp{parent.xp}
    ,yp{parent.yp}
    ,resources{0.0}
{
    std::uniform_real_distribution<double> uniform{0.0,1.0};
    std::normal_distribution<double> mutational_effect{0.0,params.sdmu};

    if (uniform(rng_r) < params.mu_x)
    {
        x += mutational_effect(rng_r);

        x = std::clamp(x,0.0,1.0);
    }
    
    if (uniform(rng_r) < params.mu_y)
    {
        y += mutational_effect(rng_r);
        y = std::clamp(x,0.0,1.0);
    }
    
    if (uniform(rng_r) < params.mu_xp)
    {
        xp += mutational_effect(rng_r);
        xp = std::clamp(xp,-10.0,10.0);
    }
    
    if (uniform(rng_r) < params.mu_yp)
    {
        yp += mutational_effect(rng_r);
        yp = std::clamp(yp,-10.0,10.0);
    }

    // give some initial bit of resources
    resources = params.resource_variation * uniform(rng_r);
}

void Individual::operator=(Individual const &other)
{
    x = other.x;
    xp = other.xp;
    y = other.y;
    yp = other.yp;

    resources = other.resources;
}
