#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>
#include "simulation.hpp"
#include "individual.hpp"
#include "parameters.hpp"


// constructor to initialize the simulation
Simulation::Simulation(Parameters const &params) :
    data_file{params.file_base_name.c_str()}
    ,rd{} // initialize random device, see *.h file
    ,seed{rd()} // initialize seed
    ,rng_r{seed} // initialize the random number generator
    ,uniform{0.0,1.0} // initialize the uniform distribution
    ,uniform_pm{-1.0,1.0} // initialize the uniform distribution +1 - 1
    ,params{params} // store parameter object
    ,singles{params.N,Individual(params.init_x,params.init_y)} // initialize pop
{} // end constructor

// run the actual simulation over max_time time steps
void Simulation::run()
{
    write_data_headers();

    for (time_step = 0; time_step <= params.max_time; ++time_step)
    {
        pair_up();

        interact();

        reproduce();

        dismiss_partner();

        mortality();


        if (time_step % params.data_interval == 0)
        {
            write_data();
        }
    }

    write_parameters();
} // end run()

// pair up single individuals
void Simulation::pair_up()
{
    // all singles pair up
    // we just do this by shuffling the list and then individuals
    // in indices 0, 1 interact. 
    //
    // final index in case of unevenness has bad luck and does not interact 
    std::shuffle(singles.begin(), singles.end(),rng_r);
} // pair_up()

// then calculate payoff according to a pd scenario
double Simulation::payoff_pd(double const x, double const xprime)
{
    return((x + xprime)/(1.0 + x + xprime) - 0.8 * (x + x * x));
}

// have individuals interact and calculate payoffs
void Simulation::interact()
{
    double x1, x2;

    // first have new pairs interact
    // loop over each pair, simply by loooping over two singles at a time
    for (int new_pair_idx = 0; 
            new_pair_idx < singles.size(); new_pair_idx += 2)
    {
        x1 = singles[new_pair_idx].x;
        x2 = singles[new_pair_idx + 1].x;

        singles[new_pair_idx].resources += 
            payoff_pd(x1,x2) - params.startup_cost;

        singles[new_pair_idx + 1].resources += 
            payoff_pd(x2,x1) - params.startup_cost;
        
    }

    for (int pair_idx = 0; pair_idx < paired.size(); pair_idx += 2)
    {
        x1 = paired[pair_idx].x;
        x2 = paired[pair_idx + 1].x;
        
        paired[pair_idx].resources += payoff_pd(x1,x2);
        paired[pair_idx + 1].resources += payoff_pd(x2,x1);
    }
}// Simulation::interact

// calculate # deaths
int Simulation::calculate_mortalities()
{
    int n = singles.size() + paired.size();

    std::binomial_distribution<int> deaths(n, params.mortality_prob);

    return(deaths(rng_r));
} // end calculate_mortalities

// produce offspring dependent on payoffs
void Simulation::reproduce()
{
    // remove offspring from previous time step
    offspring.clear();

    // make distribution of resources of all singles and paireds
    std::vector <double> resources;

    // make cumulative distributions of both singles and paired
    for (int new_pair_idx = 0; new_pair_idx < singles.size(); ++new_pair_idx)
    {
        resources.push_back(singles[new_pair_idx].resources);
    }
    
    for (int pair_idx = 0; pair_idx < paired.size(); ++pair_idx)
    {
        resources.push_back(paired[pair_idx].resources);
    }

    std::discrete_distribution<int> resource_distribution(resources.begin(), resources.end());

    // calculate mortalities
    number_mortalities = calculate_mortalities();

    int parent_idx;

    for (int offspring_idx = 0; 
            offspring_idx < number_mortalities; ++offspring_idx)
    {
        // sample parent from resource distribution
        parent_idx = resource_distribution(rng_r);

        // TODO potential for bugs here
        if (parent_idx >= singles.size())
        {
            parent_idx -= singles.size();

            assert(parent_idx >= 0);
            assert(parent_idx < paired.size());
        
            offspring.push_back(
                    Individual(paired[parent_idx]
                        ,params
                        ,rng_r));
        }
        else
        {
            assert(parent_idx >= 0);
            assert(parent_idx < singles.size());

            offspring.push_back(
                    Individual(singles[parent_idx]
                        ,params
                        ,rng_r));
        }
    } // end for offspring idx
} // end Simulation::reproduce()

void Simulation::dismiss_partner()
{
    double x1, x2;
    double y1, y2;

    std::vector <Individual> new_old_pairs;

    assert(!singles.empty());

    for (std::vector<Individual>::iterator pair_iter = singles.begin();
            pair_iter != singles.end();

        )
    {
        auto pair_iter2 = std::next(pair_iter,1);

        x1 = pair_iter->x;
        x2 = pair_iter2->x;

        y1 = pair_iter->y;
        y2 = pair_iter2->y;

        // see whether individuals will remain in pair
        if (y1 <= x2 + params.dismiss_error * uniform_pm(rng_r)
                && 
                y2 <= x1 + params.dismiss_error * uniform_pm(rng_r))
        {
            paired.push_back(*pair_iter);
            paired.push_back(*pair_iter2);


            // erase both individuals from the singles
            pair_iter = singles.erase(pair_iter);
            pair_iter = singles.erase(pair_iter);
        
        }
        else
        {
            // advance the iterator in twos (coz pairs)
            pair_iter = std::next(pair_iter,2);
        }
    }
}

// mortality events as many as there are offspring
void Simulation::mortality()
{
    // aux variab
    double prob_single;

    for (int offspring_idx = 0; offspring_idx < offspring.size(); ++offspring_idx)
    {
        prob_single = (double) singles.size() / (singles.size() + paired.size());

        if (uniform(rng_r) < prob_single)
        {
            std::uniform_int_distribution <int> single_sampler(0,singles.size() - 1);

            singles.erase(singles.begin() + single_sampler(rng_r));
        }
        else
        {
            std::uniform_int_distribution <int> paired_sampler(0,paired.size() - 1);

            int random_paired_idx = paired_sampler(rng_r);
            
            // get index of partner
            int random_paired_partner_idx = 
                random_paired_idx % 2 == 0 ? 
                    random_paired_idx + 1
                    :
                    random_paired_idx - 1;

            assert(random_paired_partner_idx >= 0);

            assert(random_paired_partner_idx < paired.size());

            // move remaining individual from this pair to stack of singles
            singles.push_back(paired[random_paired_idx + 1]);
            
            paired.erase(paired.begin() + random_paired_idx);
        }
    }
} // end Simulation::mortality()
  
void Simulation::write_data()
{
    double mean_x = 0.0;
    double ss_x = 0.0;
    double mean_y = 0.0;
    double ss_y = 0.0;

    double x,y;

    // calculate stats
    for (std::vector<Individual>::iterator single_iter = singles.begin();
            single_iter != singles.end();
            ++single_iter)
    {
        x = single_iter->x;
        mean_x += x;
        ss_x += x*x;

        y = single_iter->y;
        mean_y += y;
        ss_y += y*y;
    }

    for (std::vector<Individual>::iterator paired_iter = paired.begin();
            paired_iter != paired.end();
            ++paired_iter)
    {
        x = paired_iter->x;
        mean_x += x;
        ss_x += x*x;

        y = paired_iter->y;
        mean_y += y;
        ss_y += y*y;
    }

    int n = singles.size() + paired.size();

    mean_x /= n;
    double var_x = ss_x / n - mean_x * mean_x;

    mean_y /= n;
    double var_y = ss_y / n - mean_y * mean_y;

    data_file << time_step << ";"
        << mean_x << ";"
        << mean_y << ";"
        << var_x << ";"
        << var_y << ";"
        << paired.size() << ";"
        << singles.size() << ";"
        << std::endl;
} // end write_data();

void Simulation::write_data_headers()
{
    data_file << "time;x;y;var_x;var_y;npaired;nsingle;" << std::endl;
}

void Simulation::write_parameters()
{
    data_file << std::endl
        << std::endl
        << "mu_x;" << params.mu_x << std::endl
        << "mu_y;" << params.mu_y << std::endl
        << "mu_xp;" << params.mu_xp << std::endl
        << "mu_yp;" << params.mu_yp << std::endl
        << "sdmu;" << params.sdmu << std::endl
        << "N;" << params.N << std::endl
        << "init_x;" << params.init_x << std::endl
        << "init_y;" << params.init_y << std::endl
        << "startup_cost;" << params.startup_cost << std::endl
        << "mortality_prob;" << params.mortality_prob << std::endl
        << "dismiss_error;" << params.dismiss_error << std::endl;
}
