#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

class Parameters
{
    public:
        // probability that the 
        // interaction will be continued
        double w = 0.5;

        double mu_x = 0.01;
        double mu_y = 0.01;
        double mu_xp = 0.01;
        double mu_yp = 0.01;

        double sdmu = 0.02;
};

#endif
