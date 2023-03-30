import os
import copy
import numpy as np

# benefit function in a CPD 
def BPD(x, xother):
    return(5 * xother)

# benefit function under snowdrift
def Bsnowdrift(x, xother):
    return((x + xother) / (1 + x + xother))

def Csnowdrift(x):
    return(0.8 * (x + x*x))

class McNamara2008:

    def __init__(self, G = 30, 
                 snowdrift=False, 
                 muCoop = 0.01, 
                 muChoice = 0.01,
                 initvals = (0,0),
                 xmax = 1.0,
                 max_time = 1):

        self.G = G
        self.xmax = xmax
        self.init_vectors(initial_values = initvals)

        self.max_time = max_time

        # dynamically assign benefit and cost function
        # for this simulation dependent on what version of the model
        # we actually are interested in
        self.B = BPD;


    def run(self):

        for i in range(0, self.max_time):
            self.pairing()


    # just random pairing, choice will be done during
    # an initial round
    def pairing(self):

        sum_u = self.u.sum();
        
        # calculate mean levels in single individuals
        for i in range(0,self.G + 1):
            for j in range(0,self.G + 1):
                for k in range(0,self.G + 1):
                    for l in range(0,self.G + 1):
                        self.Q[i,j,k,l] = self.u[i,j] * self.u[k,l] / sum_u

    def init_vectors(self, initial_values):

        # initialize frequency 
        # of single individuals, u
        # rows is effort level
        # cols is choice level
        self.u = np.zeros((self.G + 1, self.G + 1))

        # calculate the indices where we need to set the frequency to 1
        # this will be the initial population

        # if trait value is given by trait = (j / G) * xmax
        # then index is given by trait * G / xmax
        index_init_effort = round(initial_values[0] * self.G / self.xmax)
        index_init_choice = round(initial_values[1] * self.G / self.xmax)

        self.u[index_init_effort,index_init_choice] = 1.0

        # initialize frequency of surviving pairs
        self.P = np.zeros((self.G + 1, self.G + 1, self.G + 1, self.G + 1))
        
        # initialize frequency of new pairs, all at 0
        self.Q = copy.deepcopy(self.P)

    def means(self):
    
        means = [0.0,0.0]

        # calculate mean levels in single individuals
        for i in range(0,self.G + 1):
            for j in range(0,self.G + 1):
                
                means[0] += self.u[i,j] * i / self.G * self.xmax
                means[1] += self.u[i,j] * j / self.G * self.xmax

                for k in range(0,self.G + 1):
                    for l in range(0,self.G + 1):
                        means[0] += (self.P[i,j,k,l] * i + self.P[i,j,k,l] * k) / self.G * self.xmax
                        means[1] += (self.P[i,j,k,l] * j + self.P[i,j,k,l] * l) / self.G * self.xmax

                        means[0] += (self.Q[i,j,k,l] * i + self.Q[i,j,k,l] * k) / self.G * self.xmax
                        means[1] += (self.Q[i,j,k,l] * j + self.Q[i,j,k,l] * l) / self.G * self.xmax

        return tuple(means)


    

