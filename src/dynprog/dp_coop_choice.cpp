#include <iostream>
#include <fstream>
#include <vector>
#include <random>

// setting up the random number 
// generator machinery
std::random_device rd;
unsigned seed = rd();
std::mt19937 mt(seed); // random number generator

// then prepare the (0,1)-uniform distribution
// which we need to translate probablities into
// actual events happening.
std::uniform_real_distribution<double> uniform(0, 1); // real number between 0 and 1 (uniform)

double init_threshold{0.0};

std::string file_name{"data_coop_choice_"};

// max resources
int max_resources{100};

// threshold as a function 
// of resource dependence
// this is a vector for each level
// of resources
std::vector <double> threshold;

// investment as a function
// of resources
// this
std::vector <double> investment;

// obtain some parameter values from the command
// line so that we do not have to recompile the programm
// every time we want to run it for a different set of
// parameters
void get_command_line_args(int argc, char *argv[])
{
    max_resources = std::stoi(argv[1]);
    file_name = argv[2];
} // end command_line_args()

// initialize all variables
void initialize_variables()
{
    // make a dummy threshold and assign it a uniform
    // value 
    std::vector <double> init_threshold_vector(max_resources, init_threshold);

    threshold = init_threshold_vector;

    std::vector <double> init_investment_vector(max_resources, init_investment);

    investment = init_investment_vector;

} // end initialize_variables()

int main(int argc, char *argv[])
{
    get_command_line_args(argc, argv);
    initialize_variables();
} // end int main()
