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
    initial_values = [0.3,0.3] # first elmt coop, next one choice
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


    def run(self):

        for i in range(0, self.params.max_time):
            self.pairing()
            self.payoffs()
            self.reproduce()
            self.dismissal()

    def B(self, i, j):
        return((i + j)/(1 + i + j))

    def C(self, i):
        return(0.8 * (i + i*i))

    def W(self, i, j):

        return(self.B(i,j) - self.C(i))

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

        sum_u = self.u.sum();
        
        # calculate mean levels in single individuals
        for i in range(0,self.params.G + 1):
            for j in range(0,self.params.G + 1):
                for k in range(0,self.params.G + 1):
                    for l in range(0,self.params.G + 1):
                        self.Q[i,j,k,l] = self.u[i,j] * self.u[k,l] / sum_u

    # dismiss incompatible partners
    def dismissal(self):
        
        print("Q pre dismissal.")
        print(self.Q.sum())
        
        print("P pre dismissal.")
        print(self.P.sum())

        for i in range(0,self.params.G + 1):
            for j in range(0,self.params.G + 1):
                for k in range(0,self.params.G + 1):
                    for l in range(0,self.params.G + 1):
                        if (i >= l and k >= j):
                            self.P[i,j,k,l] += self.Q[i,j,k,l]

                            self.Q[i,j,k,l] = 0.0
                            
        print("Q post dismissal.")
        print(self.Q.sum())
        
        print("P post dismissal.")
        print(self.P.sum())

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

        # now locate this somewhere in a grid of evenly spaced values of 0 - G
        histogram_effort = np.histogram(
                a=[initial_effort_phenotype],
                bins=self.params.G + 1,
                range=(0,self.params.xmax))[0]

        print(histogram_effort)
        sys.exit(1)

        print(index_init_effort)
        print(index_init_choice)

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


    

