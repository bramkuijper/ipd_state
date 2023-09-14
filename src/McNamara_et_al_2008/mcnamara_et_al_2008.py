import os
import copy
import sys
import numpy as np
from enum import IntEnum

# benefit function in a continuous prisoner's dilemma
def BPD(x, xother):
    return(5 * xother)

# benefit function under snowdrift
def Bsnowdrift(x, xother):
    return((x + xother) / (1 + x + xother))

def Csnowdrift(x):
    return(0.8 * (x + x*x))

# Enum type class to distinguish between both traits
class Traits(IntEnum):
    Investment = 0
    Choice = 1


class Parameters:

    G = 30
    snowdrift = False

    # values between 0 and 1
    initial_values = [0.05,0.05] # first elmt coop, next one choice
    max_time = 500
    xmax = 0.18
    
    mu = [0.05,0.05] # first elmt coop, next one choice

    A = 0.01
    S = 0.01

    M = 0.01


# set up the class corresponding to this simulation
class McNamara2008:

    # constructor with default param values
    def __init__(self, parameters):

        self.params = parameters
        self.n_alleles = self.params.G + 1

        self.init_vectors()

    def run(self):

        self.output_data_headers()

        self.fillWmatrix()

        for self.time_step in range(0, self.params.max_time):
            self.pairing()
            self.payoffs()
            self.reproduce()
            self.dismissal_mortality()
            self.output_data()

        self.output_parameters()

    def B(self, x_i, x_j):
        return((x_i + x_j)/(1 + x_i + x_j))

    def C(self, x_i):
        return(0.8 * (x_i + x_i*x_i))

    def W(self, i, j):
        
        x_i = i/self.n_alleles * self.params.xmax
        x_j = j/self.n_alleles * self.params.xmax
        return(self.B(x_i,x_j) - self.C(x_i))

    def fillWmatrix(self):
        
        self.Wval = np.zeros((self.n_alleles,self.n_alleles))
        
        for i in range(0,self.n_alleles):
            for j in range(0,self.n_alleles):

                self.Wval[i,j] = self.W(i,j)



    def mutate(self, k, i, trait):

        if abs(k - i) > 1:
            return(0)

        if i == 0: 

            # return 1-mu/2 if i == k
            if k == 0:
                return(1.0 - self.params.mu[trait]/2.0)
            
            # otherwise k necessarily 1, return mu/2
            return(self.params.mu[trait]/2.0)

        if i == k:
            return(1.0 - self.params.mu[trait]);

        return(self.params.mu[trait]/2.0)


    # pairing of single individual
    # just random pairing, choice will be done during
    # an initial round
    def pairing(self):

        self.sum_u = self.u.sum();


        self.Q = np.zeros(
                (
                    self.n_alleles,
                    self.n_alleles,
                    self.n_alleles,
                    self.n_alleles
                    ))
        
        # calculate mean levels in single individuals
        for i in range(0,self.n_alleles):
            for j in range(0,self.n_alleles):
                for k in range(0,self.n_alleles):
                    for l in range(0,self.n_alleles):
                        self.Q[i,j,k,l] = self.u[i,j] * self.u[k,l] / self.sum_u

    # calculate rho in eq. (A1)
    def calculate_rho(self):

        for i in range(0,self.n_alleles):
            for j in range(0,self.n_alleles):
                
                sum_P = 0.0
                
                for k in range(0,self.n_alleles):
                    for l in range(0,self.n_alleles):

                        sum_P += self.P[i,j,k,l]

                self.rho[i,j] = sum_P + self.u[i,j]


    # dismiss incompatible partners and have some semblance of mortality
    def dismissal_mortality(self):
        
        # survival probability of a pair
        prob_pair_survival = (1.0 - self.params.M) * (1.0 - self.params.M)

        # calculate rho
        self.calculate_rho()



        a = self.uprime = np.zeros((self.n_alleles, self.n_alleles))

        for i in range(0,self.n_alleles):
            for j in range(0,self.n_alleles):

                sumPprime = 0.0

                for k in range(0,self.n_alleles):
                    for l in range(0,self.n_alleles):

                        self.P[i,j,k,l] = prob_pair_survival * self.P[i,j,k,l]

                        # no dismissal
                        if (i >= l and k >= j):
                            self.P[i,j,k,l] += prob_pair_survival * self.Q[i,j,k,l]

                        sumPprime += self.P[i,j,k,l]

                # equation (A7)
                a[i,j] = (1.0 - self.params.M) * self.rho[i,j] - sumPprime 

                # equation (A8)
                self.uprime[i,j] = a[i,j] + self.params.M * self.v[i,j] / self.V
                

        # now normalize u and P
        self.u = self.uprime / self.uprime.sum()
        self.P = self.P / self.P.sum()

        # we have already updated P to P' in the loop
        assert round(self.P.sum(),5) >= 0
        assert round(self.P.sum(),5) <= 1
        assert round(self.u.sum(),5) >= 0
        assert round(self.u.sum(),5) <= 1


    def output_data_headers(self):

        header_str = "time;"

        for i in range(0,self.n_alleles):
            header_str += f"i{i};c{i};"

        header_str += "mean_i;mean_c;V;"

        print(header_str)


    def output_data(self):

        mean_x = 0.0
        mean_y = 0.0

        xmax_per_allele = 1.0 / self.n_alleles  * self.params.xmax

        # marginal frequencies
        i_marg = np.zeros((self.n_alleles))
        j_marg = np.zeros((self.n_alleles))

        for i in range(0, self.n_alleles):
            for j in range(0, self.n_alleles):
                mean_x += self.rho[i,j] * i * xmax_per_allele
                mean_y += self.rho[i,j] * j * xmax_per_allele

                i_marg[i] += self.rho[i,j]
                j_marg[j] += self.rho[i,j]

        # then print the stuff
        outputstr = f"{self.time_step};"

        for i in range(0, self.n_alleles):
            outputstr += f"{i_marg[i]};{j_marg[i]};"

        outputstr += f"{mean_x};{mean_y};{self.V};"

        print(outputstr)

    # output a two-dimensional matrix for error checking purposes
    def output_matrix(self, the_matrix):

        dimensions = the_matrix.shape

        output_str = ""

        for i in range(0, dimensions[0]):

            if i == 0:
                # print columns

                # empty entry coz top-left 
                # corner of matrix
                output_str += ";"

                for j in range(0, dimensions[1]):
                    output_str += f"{j};"

                output_str += "\n"

            output_str += f"{i};"

            for j in range(0, dimensions[1]):
                output_str += f"{the_matrix[i,j]}" + ";"

            output_str += "\n"

        print(output_str)

    def output_parameters(self):
        output_str ="\n\n"

        output_str += f"G;{self.params.G}\n"

        output_str += f"init_i;{self.params.initial_values[0]}\n"
        output_str += f"init_c;{self.params.initial_values[1]}\n"
        
        output_str += f"xmax;{self.params.xmax}\n"
        
        output_str += f"mu_i;{self.params.mu[0]}\n"
        output_str += f"mu_c;{self.params.mu[1]}\n"
    
        output_str += f"A;{self.params.A}\n"
        output_str += f"S;{self.params.S}\n"
        output_str += f"M;{self.params.M}\n"

        print(output_str)


    # calculate payoffs
    def payoffs(self):

        # calculate mean levels in single individuals
        for i in range(0,self.n_alleles):
            for j in range(0,self.n_alleles):

                # reset payoffs 
                self.r[i,j] = 0.0

                for k in range(0,self.n_alleles):
                    for l in range(0,self.n_alleles):

                        self.r[i,j] += self.P[i,j,k,l] * (self.Wval[i,k] + self.params.A) + self.Q[i,j,k,l] * (self.Wval[i,k] + self.params.A - self.params.S)

    # reproduce
    def reproduce(self):


        for k in range(0,self.n_alleles):
            for l in range(0,self.n_alleles):
                self.v[k,l] = 0.0

                for i in range(0,self.n_alleles):
                    for j in range(0,self.n_alleles):
                        self.v[k,l] += self.r[i,j] * self.mutate(k, i, Traits.Investment) * self.mutate(l, j, Traits.Choice)


        # calculate the full sum
        self.V = self.v.sum()

    def phenotype_to_genotype(self, phenotype):
        # now locate this phenotype in a grid of evenly spaced values of 0 - G
        the_hist = list(np.histogram(
                a=[phenotype],
                bins=self.n_alleles,
                range=(0,self.params.xmax))[0])

        # this should return a list of 0s and a single 1
        # return the index of that 1
        return(the_hist.index(1))


    # initialize all the vectors
    def init_vectors(self):

        # initialize frequency 
        # of single individuals, u
        # rows is effort level
        # cols is choice level
        self.u = self.rho = np.zeros((self.n_alleles, self.n_alleles))

        # calculate the indices where we need to set the frequency to 1
        # this will be the initial population
        initial_effort_phenotype = self.params.initial_values[0] * self.params.xmax
        initial_choice_phenotype = self.params.initial_values[1] * self.params.xmax

        # from that initial effort phenotype, calculate the gene index
        index_init_effort = self.phenotype_to_genotype(initial_effort_phenotype)
        index_init_choice = self.phenotype_to_genotype(initial_choice_phenotype)

        self.u[index_init_effort,index_init_choice] = 1.0

        # initialize frequency of surviving pairs
        self.P = np.zeros((self.n_alleles, self.n_alleles, self.n_alleles, self.n_alleles))
        
        # initialize frequency of new pairs, all at 0
        self.Q = copy.deepcopy(self.P)

        # vector for payoffs
        self.r = np.zeros((self.n_alleles, self.n_alleles))

        # vector for payoffs x mutation rates
        self.v = np.zeros((self.n_alleles, self.n_alleles))

        # total V
        self.V = 0.0;

    def means(self):
    
        means = [0.0,0.0]

        # calculate mean levels of choice and cooperation 
        # in single individuals
        for i in range(0,self.n_alleles):
            for j in range(0,self.n_alleles):
                
                means[0] += self.u[i,j] * i / self.params.G * self.params.xmax
                means[1] += self.u[i,j] * j / self.params.G * self.params.xmax

                for k in range(0,self.n_alleles):
                    for l in range(0,self.n_alleles):
                        means[0] += (self.P[i,j,k,l] * i + self.P[i,j,k,l] * k) / self.params.G * self.params.xmax
                        means[1] += (self.P[i,j,k,l] * j + self.P[i,j,k,l] * l) / self.params.G * self.params.xmax

                        means[0] += (self.Q[i,j,k,l] * i + self.Q[i,j,k,l] * k) / self.params.G * self.params.xmax
                        means[1] += (self.Q[i,j,k,l] * j + self.Q[i,j,k,l] * l) / self.params.G * self.params.xmax

        return tuple(means)


    

