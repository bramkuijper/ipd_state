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
    mu = [0.01,0.01] # first elmt coop, next one choice

    # values between 0 and 1
    initial_values = [0.3,0.1] # first elmt coop, next one choice
    max_time = 1
    xmax = 0.18

    A = 0.01
    S = 0.01


# set up the class corresponding to this simulation
class McNamara2008:

    # constructor with default param values
    def __init__(self, parameters):

        self.params = parameters
        self.init_vectors()

        self.n_alleles = self.params.G + 1


    def run(self):

        for i in range(0, self.params.max_time):
            self.pairing()
            self.payoffs()
            self.reproduce()
            self.dismissal()

    def B(self, x_i, x_j):
        return((x_i + x_j)/(1 + x_i + x_j))

    def C(self, x_i):
        return(0.8 * (x_i + x_i*x_i))

    def W(self, i, j):
        
        x_i = i/self.n_alleles * self.params.xmax
        x_j = j/self.n_alleles * self.params.xmax
        return(self.B(x_i,x_j) - self.C(x_i))

    def mutate(self, k, i, trait):

        if abs(k - i) > 1:
            return(0)

        if i == 0: 

            # return 1-mu/2 if i == k
            if k == 0:
                return(1.0 - self.params.mu[trait]/2)
            
            # otherwise k necessarily 1, return mu/2
            return(self.params.mu[trait]/2)

        if i == k:
            return(1.0 - self.params.mu[trait]);

        return(self.params.mu[trait]/2)


    # pairing of single individual
    # just random pairing, choice will be done during
    # an initial round
    def pairing(self):

        self.sum_u = self.u.sum();
        
        # calculate mean levels in single individuals
        for i in range(0,self.params.G + 1):
            for j in range(0,self.params.G + 1):
                for k in range(0,self.params.G + 1):
                    for l in range(0,self.params.G + 1):
                        self.Q[i,j,k,l] = self.u[i,j] * self.u[k,l] / self.sum_u

        print("sum u:")
        print(self.sum_u)

    # dismiss incompatible partners
    def dismissal(self):
        
        print("Q pre dismissal.")
        print(self.Q.sum())
        
        print("P pre dismissal.")
        print(self.P.sum())

        # store the original vector with u frequencies
        # we use this to calculate the new u frequencies
        # from Q
        u_tmp = copy.deepcopy(self.u)
        sum_u_tmp = self.sum_u

        self.u = np.zeros((self.params.G + 1, self.params.G + 1))

        for i in range(0,self.params.G + 1):
            for j in range(0,self.params.G + 1):

                for k in range(0,self.params.G + 1):
                    for l in range(0,self.params.G + 1):

                        # no dismissal
                        if (i >= l and k >= j):
                            self.P[i,j,k,l] += self.Q[i,j,k,l]
                        else: # dismissal, add values to the u vector agin
                            print(f"{i} {j} {k} {l} {u_tmp[k,l]} {self.Q[i,j,k,l]}")
                            self.u[i,j] += self.Q[i,j,k,l] / u_tmp[k,l] * sum_u_tmp;
                            self.u[k,l] += self.Q[i,j,k,l] / u_tmp[i,j] * sum_u_tmp;
        
        print("Q (pre-paired) post dismissal.")
        print(self.Q.sum())
        
        print("P (paired) post dismissal.")
        print(self.P.sum())

        print("u (unpaired) post dismissal.")
        print(self.u.sum())

        sys.exit(1)


    # calculate payoffs
    def payoffs(self):

        # calculate mean levels in single individuals
        for i in range(0,self.params.G + 1):
            for j in range(0,self.params.G + 1):

                # reset payoffs 
                self.r[i,j] = 0.0

                for k in range(0,self.params.G + 1):
                    for l in range(0,self.params.G + 1):
                        self.r[i,j] += self.P[i,j,k,l] * (self.W(i,k) + self.params.A) + self.Q[i,j,k,l] * (self.W(i,k) + self.params.A - self.params.S)

    # reproduce
    def reproduce(self):

        for k in range(0,self.params.G + 1):
            for l in range(0,self.params.G + 1):
                self.v[k,l] = 0.0

                for i in range(0,self.params.G + 1):
                    for j in range(0,self.params.G + 1):
                        self.v[k,l] += self.r[i,j] * self.mutate(k, i, Traits.Investment) * self.mutate(l, j, Traits.Choice)

        # calculate the full sum
        self.V = self.v.sum()

    def phenotype_to_genotype(self, phenotype):
        # now locate this phenotype in a grid of evenly spaced values of 0 - G
        the_hist = list(np.histogram(
                a=[phenotype],
                bins=self.params.G + 1,
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
        self.u = np.zeros((self.params.G + 1, self.params.G + 1))

        # calculate the indices where we need to set the frequency to 1
        # this will be the initial population
        initial_effort_phenotype = self.params.initial_values[0] * self.params.xmax
        initial_choice_phenotype = self.params.initial_values[1] * self.params.xmax

        # from that initial effort phenotype, calculate the gene index
        index_init_effort = self.phenotype_to_genotype(initial_effort_phenotype)
        index_init_choice = self.phenotype_to_genotype(initial_choice_phenotype)

        self.u[index_init_effort,index_init_choice] = 1.0

        # initialize frequency of surviving pairs
        self.P = np.zeros((self.params.G + 1, self.params.G + 1, self.params.G + 1, self.params.G + 1))
        
        # initialize frequency of new pairs, all at 0
        self.Q = copy.deepcopy(self.P)

        # vector for payoffs
        self.r = np.zeros((self.params.G + 1, self.params.G + 1))

        # vector for payoffs x mutation rates
        self.v = np.zeros((self.params.G + 1, self.params.G + 1, self.params.G + 1, self.params.G + 1))

        # total V
        self.V = 0.0;

    def means(self):
    
        means = [0.0,0.0]

        # calculate mean levels of choice and cooperation 
        # in single individuals
        for i in range(0,self.params.G + 1):
            for j in range(0,self.params.G + 1):
                
                means[0] += self.u[i,j] * i / self.params.G * self.params.xmax
                means[1] += self.u[i,j] * j / self.params.G * self.params.xmax

                for k in range(0,self.params.G + 1):
                    for l in range(0,self.params.G + 1):
                        means[0] += (self.P[i,j,k,l] * i + self.P[i,j,k,l] * k) / self.params.G * self.params.xmax
                        means[1] += (self.P[i,j,k,l] * j + self.P[i,j,k,l] * l) / self.params.G * self.params.xmax

                        means[0] += (self.Q[i,j,k,l] * i + self.Q[i,j,k,l] * k) / self.params.G * self.params.xmax
                        means[1] += (self.Q[i,j,k,l] * j + self.Q[i,j,k,l] * l) / self.params.G * self.params.xmax

        return tuple(means)


    

