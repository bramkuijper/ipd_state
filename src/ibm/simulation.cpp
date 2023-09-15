#include <iostream>
#include <vector>
#include <random>
#include <cmath>
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
{
} // end constructor

// run the actual simulation over max_time time steps
void Simulation::run()
{
    write_parameters();
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

} // end run()

// pair up single individuals
void Simulation::pair_up()
{
    check_state();
    // all singles pair up
    // we just do this by shuffling the list and then individuals
    // in indices 0, 1 interact. 
    //
    // final index in case of unevenness has bad luck and does not interact 
    std::shuffle(singles.begin(), singles.end(),rng_r);

    check_state();
} // pair_up()

// then calculate payoff according to a pd scenario
double Simulation::payoff_pd(double const x, double const xprime)
{
    double xsum = x + xprime;
    double B = -1.4 * xsum * xsum + 6 * xsum;

    double C = -1.6 * x * x + 4.56 * x; 

    assert(std::isnormal(B-C));
    return(B - C);
}

double Simulation::payoff_snowdrift(double const x, double const xprime)
{
    double retval = (x + xprime)/(1.0 + x + xprime) - 0.8 * (x + x * x);

    if (retval != 0.0)
    {
        assert(std::isnormal(retval));
    }

    return(retval);
}

double Simulation::payoff(double const x, double const xprime)
{
    if (params.cooperation_type == prisoners_dilemma)
    {
        return(payoff_pd(x, xprime));
    }

    return(payoff_snowdrift(x, xprime));
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
        check_state();
        
        if (new_pair_idx + 1 >= singles.size())
        {
            break;
        }

        assert(std::fabs(singles[new_pair_idx].resources) <= params.max_resources);
        assert(std::fabs(singles[new_pair_idx + 1].resources) <= params.max_resources);

        // gift by individual 1: baseline (x) + xp * (individual 1's resources)
        x1 = singles[new_pair_idx].x + singles[new_pair_idx].xp * singles[new_pair_idx].resources;
        // gift by individual 2: baseline (x) + xp * (individual 2's resources)
        x2 = singles[new_pair_idx + 1].x + singles[new_pair_idx + 1].xp * singles[new_pair_idx + 1].resources;

        // now calculate the payoff according to PD logic
        singles[new_pair_idx].resources += 
            payoff(x1,x2) - params.startup_cost;

        singles[new_pair_idx].resources = std::clamp(
                singles[new_pair_idx].resources,0.0,params.max_resources);

        if (singles[new_pair_idx].resources > params.max_resources)
        {
            singles[new_pair_idx].resources = params.max_resources;
        }

        singles[new_pair_idx + 1].resources += 
            payoff(x2,x1) - params.startup_cost;
        
        singles[new_pair_idx + 1].resources = std::clamp(
                singles[new_pair_idx + 1].resources,0.0,params.max_resources);
        
        if (singles[new_pair_idx + 1].resources > params.max_resources)
        {
            singles[new_pair_idx + 1].resources = params.max_resources;
        }
    }

    assert(paired.size() % 2 == 0);

    for (int pair_idx = 0; pair_idx < paired.size(); pair_idx += 2)
    {
        check_state();
        assert(std::fabs(paired[pair_idx].resources) <= params.max_resources);
        assert(std::fabs(paired[pair_idx + 1].resources) <= params.max_resources);

        x1 = paired[pair_idx].x + 
            paired[pair_idx].xp * paired[pair_idx].resources;

        x2 = paired[pair_idx + 1].x + 
            paired[pair_idx + 1].xp * paired[pair_idx + 1].resources;
        
        paired[pair_idx].resources += payoff(x1,x2);
        
        paired[pair_idx].resources = std::clamp(
                paired[pair_idx].resources,0.0,params.max_resources);

        if (paired[pair_idx].resources > params.max_resources)
        {
            paired[pair_idx].resources = params.max_resources;
        }

        paired[pair_idx + 1].resources += payoff(x2,x1);
        
        paired[pair_idx + 1].resources = std::clamp(
                paired[pair_idx + 1].resources,0.0,params.max_resources);
        
        if (paired[pair_idx + 1].resources > params.max_resources)
        {
            paired[pair_idx + 1].resources = params.max_resources;
        }
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
    check_state();
    // remove offspring from previous time step
    offspring.clear();

    check_state();
    // make distribution of resources of all singles and paireds
    std::vector <double> resources;

    // make cumulative distributions of both singles and paired
    for (int new_pair_idx = 0; new_pair_idx < singles.size(); ++new_pair_idx)
    {
        assert(std::fabs(singles[new_pair_idx].resources) <= params.max_resources);
        resources.push_back(singles[new_pair_idx].resources);
    }
    
    for (int pair_idx = 0; pair_idx < paired.size(); ++pair_idx)
    {
//        std::cout << time_step << " " << pair_idx << " " << paired.size() << " " << paired[pair_idx].resources << std::endl;
        assert(std::fabs(paired[pair_idx].resources) <= params.max_resources);
        resources.push_back(paired[pair_idx].resources);
    }
    
    check_state();

    // make a distribution of all individuals 
    // according to their resources
    std::discrete_distribution<int> resource_distribution(
            resources.begin(), resources.end());

    // calculate mortalities
    number_mortalities = calculate_mortalities();

    int parent_idx;

    check_state();
    for (int offspring_idx = 0; 
            offspring_idx < number_mortalities; ++offspring_idx)
    {

        // sample parent from resource distribution
        parent_idx = resource_distribution(rng_r);
       
        // parent_idx exceeds the stack of singles?
        // hence it must be paired
        if (parent_idx >= singles.size())
        {
            parent_idx -= singles.size();

            assert(parent_idx >= 0);
            assert(parent_idx < paired.size());
            
            if (std::fabs(paired[parent_idx].resources) > params.max_resources)
            {
                std::cout << "not good: " << " " << time_step << " "  
                    << parent_idx << " " 
                    << paired[parent_idx].resources << " " 
                    << std::fabs(paired[parent_idx].resources) << " " 
                    << params.max_resources << " " 
                    << singles.size() << " " 
                    << paired.size() << " "  << std::endl;
            }
            
            assert(std::fabs(paired[parent_idx].resources) <= params.max_resources);

            // make an offspring
            Individual Kid(paired[parent_idx]
                        ,params
                        ,rng_r);
        
            offspring.push_back(Kid);
        
            assert(std::fabs(paired[parent_idx].resources) <= params.max_resources);
            assert(Kid.resources <= params.max_resources);

            paired[parent_idx].resources -= params.cost_of_reproduction;

            if (paired[parent_idx].resources < 0)
            {
                paired[parent_idx].resources = 0.0;
            }
        }
        else
        {
            assert(parent_idx >= 0);
            assert(parent_idx < singles.size());

            Individual Kid(singles[parent_idx]
                        ,params
                        ,rng_r);

            offspring.push_back(Kid);

            assert(std::fabs(singles[parent_idx].resources) <= params.max_resources);
            assert(Kid.resources <= params.max_resources);
            singles[parent_idx].resources -= params.cost_of_reproduction;
            
            if (singles[parent_idx].resources < 0)
            {
                singles[parent_idx].resources = 0.0;
            }
        }
    } // end for offspring idx
    check_state();
} // end Simulation::reproduce()

void Simulation::dismiss_partner()
{
    double x1, x2;
    double y1, y2;

    std::vector <Individual> new_old_pairs;

    if (singles.empty())
    {
        return;
    }

    check_state();
//    int i = 0;

    //  loop through newly formed pairs
    for (std::vector<Individual>::iterator pair_iter = singles.begin();
            pair_iter != singles.end();

        )
    {
//        std::cout << time_step << " " << "pair " << i << " " << singles.size() << " " << std::distance(pair_iter, singles.end()) << " " << std::endl;
//        ++i;
        if (std::distance(pair_iter, singles.end()) <= 1)
        {
            break;
        }

        auto pair_iter2 = std::next(pair_iter,1);

        assert(pair_iter->resources <= params.max_resources);
        assert(pair_iter2->resources <= params.max_resources);
        x1 = pair_iter->x + pair_iter->xp * pair_iter->resources;
        x2 = pair_iter2->x + pair_iter2->xp * pair_iter2->resources;

        y1 = pair_iter->y + pair_iter->yp * pair_iter->resources;
        y2 = pair_iter2->y + pair_iter2->yp * pair_iter2->resources;

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
            if (std::distance(pair_iter, singles.end()) <= 2)
            {
                break;
            }
            // advance the iterator in twos (coz pairs)
            pair_iter = std::next(pair_iter,2);
        }
    }
    
    check_state();
} // end void Simulation::dismiss_partner

// Bounds checking on all the individual's states
void Simulation::check_state()
{
    for (std::vector <Individual>::iterator ind_iter = singles.begin();
            ind_iter != singles.end();
            ++ind_iter)
    {
        assert(std::fabs(ind_iter->resources) <= params.max_resources);
        assert(ind_iter->y <= 1.0);
    }

    for (std::vector <Individual>::iterator ind_iter = paired.begin();
            ind_iter != paired.end();
            ++ind_iter)
    {
        assert(std::fabs(ind_iter->resources) <= params.max_resources);
        assert(ind_iter->y <= 1.0);
    }

} // end check_state()

// mortality events as many as there are offspring
void Simulation::mortality()
{
    check_state();

    // aux variab
    double prob_single;

    // have offspring replace existing breeders
    for (int offspring_idx = 0; offspring_idx < offspring.size(); ++offspring_idx)
    {
        prob_single = (double) singles.size() / (singles.size() + paired.size());

        assert(prob_single >= 0.0);
        assert(prob_single <= 1.0);

        if (uniform(rng_r) < prob_single)
        {
            assert(singles.size() > 0);
            std::uniform_int_distribution <int> single_sampler(0,singles.size() - 1);

            singles.erase(singles.begin() + single_sampler(rng_r));
        }
        else
        {
            assert(paired.size() >= 2);
            assert(paired.size() % 2 == 0);

            // death of a pair member
            std::uniform_int_distribution <int> paired_sampler(0,paired.size() - 1);

            int random_paired_idx = paired_sampler(rng_r);
            
            // get index of partner
            // if even, i.e., 0, 2, 4, etc, we have the first member of a pair
            // hence the partner index is the focal + 1
            // other the partner's index is the focal's - 1
            if (random_paired_idx % 2 == 0)
            {
                assert(random_paired_idx + 1 > 0);
                assert(random_paired_idx + 1 < paired.size());

                // partner is in position random_paired_idx + 1
                singles.push_back(paired[random_paired_idx + 1]);

                // remove focal because of mortality
                paired.erase(paired.begin() + random_paired_idx);

                // remove focal's partner
                paired.erase(paired.begin() + random_paired_idx);
            }
            else 
            {
                assert(random_paired_idx - 1 >= 0);

                assert(random_paired_idx - 1 < paired.size());

                // partner is in position random_paired_idx + 1
                singles.push_back(paired[random_paired_idx - 1]);

                // remove focal because of mortality
                paired.erase(paired.begin() + random_paired_idx - 1);

                // remove focal's partner
                paired.erase(paired.begin() + random_paired_idx - 1);

            }
        }
    }

    check_state();

    // add offspring to the stack of singles
    singles.insert(singles.end(), offspring.begin(), offspring.end());

    check_state();
} // end Simulation::mortality()
  
void Simulation::write_data()
{
    double mean_x = 0.0;
    double ss_x = 0.0;
    double mean_y = 0.0;
    double ss_y = 0.0;

    double mean_xp = 0.0;
    double ss_xp = 0.0;
    double mean_yp = 0.0;
    double ss_yp = 0.0;

    double x,y;

    double mean_resources_single = 0.0;
    double mean_resources_paired = 0.0;
    double mean_resources = 0.0;
    double ss_resources_single = 0.0;
    double ss_resources_paired = 0.0;
    double ss_resources = 0.0;

    // calculate stats
    for (std::vector<Individual>::iterator single_iter = singles.begin();
            single_iter != singles.end();
            ++single_iter)
    {
        x = single_iter->x;
        mean_x += x;
        ss_x += x*x;

        x = single_iter->xp;
        mean_xp += x;
        ss_xp += x*x;

        y = single_iter->y;
        mean_y += y;
        ss_y += y*y;
        
        y = single_iter->yp;
        mean_yp += y;
        ss_yp += y*y;

        assert(single_iter->resources <= params.max_resources);
        x = single_iter->resources;
        mean_resources_single += x;
        ss_resources_single += x*x;

        mean_resources += x;
        ss_resources += x*x;
    }

    for (std::vector<Individual>::iterator paired_iter = paired.begin();
            paired_iter != paired.end();
            ++paired_iter)
    {
        x = paired_iter->x;
        mean_x += x;
        ss_x += x*x;

        x = paired_iter->xp;
        mean_xp += x;
        ss_xp += x*x;

        y = paired_iter->y;
        mean_y += y;
        ss_y += y*y;
        
        y = paired_iter->yp;
        mean_yp += y;
        ss_yp += y*y;
        
        assert(paired_iter->resources <= params.max_resources);
        x = paired_iter->resources;
        mean_resources_paired += x;
        ss_resources_paired += x*x;

        mean_resources += x;
        ss_resources += x*x;
    }

    int n = singles.size() + paired.size();

    mean_x /= n;
    double var_x = ss_x / n - mean_x * mean_x;

    mean_xp /= n;
    double var_xp = ss_xp / n - mean_xp * mean_xp;

    mean_y /= n;
    double var_y = ss_y / n - mean_y * mean_y;
    
    mean_yp /= n;
    double var_yp = ss_yp / n - mean_yp * mean_yp;

    mean_resources /= n;
    double var_resources = ss_resources / n - mean_resources * mean_resources;
    
    mean_resources_single /= singles.size();
    double var_resources_single = ss_resources_single / singles.size() - mean_resources_single * mean_resources_single;
    
    mean_resources_paired /= paired.size();
    double var_resources_paired = ss_resources_paired / paired.size() - mean_resources_paired * mean_resources_paired;

    data_file << time_step << ";"
        << mean_x << ";"
        << mean_y << ";"
        << mean_xp << ";"
        << mean_yp << ";"
        << var_x << ";"
        << var_y << ";"
        << var_xp << ";"
        << var_yp << ";"
        << mean_resources << ";"
        << var_resources << ";"
        << mean_resources_single << ";"
        << var_resources_single << ";"
        << mean_resources_paired << ";"
        << var_resources_paired << ";"
        << paired.size() << ";"
        << singles.size() << ";"
        << number_mortalities << ";"
        << std::endl;

} // end write_data();

void Simulation::write_data_headers()
{
    data_file << 
        "time;x;y;xp;yp;var_x;var_y;var_xp;var_yp;"
        << "mean_resources;" 
        << "var_resources;" 
        << "mean_resources_single;" 
        << "var_resources_single;" 
        << "mean_resources_paired;" 
        << "var_resources_paired;" 
        << "npaired;nsingle;nmort;" << std::endl;
} // end write_data_headers()

void Simulation::write_parameters()
{
    data_file << std::endl
        << std::endl
        << "mu_x;" << params.mu_x << std::endl
        << "mu_y;" << params.mu_y << std::endl
        << "mu_xp;" << params.mu_xp << std::endl
        << "mu_yp;" << params.mu_yp << std::endl
        << "sdmu;" << params.sdmu << std::endl
        << "seed;" << seed << std::endl
        << "game_type;" << params.cooperation_type << std::endl
        << "max_resources;" << params.max_resources << std::endl
        << "resource_variation;" << params.resource_variation << std::endl
        << "cost_of_reproduction;" << params.cost_of_reproduction << std::endl
        << "N;" << params.N << std::endl
        << "init_x;" << params.init_x << std::endl
        << "init_y;" << params.init_y << std::endl
        << "startup_cost;" << params.startup_cost << std::endl
        << "mortality_prob;" << params.mortality_prob << std::endl
        << "dismiss_error;" << params.dismiss_error << std::endl;
}
