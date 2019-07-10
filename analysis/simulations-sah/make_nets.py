### Import modules

import networkx as nx
import random_modular_generator_variable_modules as rmg
import sequence_generator as sg
from sklearn.model_selection import ParameterGrid
import numpy as np
import os
import csv
import random

os.chdir('/Users/paigemiller/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-sah')

##### Parameters for SIR simulations
# Network params
N= [1000]                                    # network sizes
Q= [0, 0.15, 0.3, 0.45]                    # network modularity 
d= 10                                      # mean degree
m= 2                                       # number of modules (sexes)
sfunction = sg.scalefree_sequence          # degree distribution 
modfunction = sg.regular_sequence          # module size distribution

var_grid = list(ParameterGrid({'N' : N, 'Q' : Q}))

# Loop params
reps= 1

##### Create and save networks  ########

for x in range(0, len(var_grid)): 
    for y in range(0, reps):
        n=var_grid[x]["N"]
        q=var_grid[x]["Q"]
    
        # generating network
        G = rmg.generate_modular_networks(n, sfunction, modfunction, q, m, d)
        nx.write_graphml(G, "networks/G_Q"+str(q)+"_N"+str(n)+"_rep"+str(y)+".graphml")

        
