### Import modules

import networkx as nx
import random_modular_generator_variable_modules as rmg
import sequence_generator as sg
from sklearn.model_selection import ParameterGrid
import numpy as np
import os
import csv
import random

#os.chdir('/Users/paigemiller/Documents/UGA/phd/research-projects/miller-tb-assortativity/analysis/simulations-sah')

##### Parameters for SIR simulations
# Network params
N= [1000]                                  # network sizes
Q= [.45, .3, .15, 0]                    # network modularity 
d= 10                                      # mean degree
m= 2                                       # number of modules (sexes)
sfunction = sg.geometric_sequence          # degree distribution 
modfunction = sg.regular_sequence          # module size distribution

var_grid = list(ParameterGrid({'N' : N, 'Q' : Q}))

# Loop params
reps= 300

##### Create and save networks  ########

try: 
	for x in range(0, len(var_grid)): 
	    for y in range(0, 35):
	        n=var_grid[x]["N"]
	        q=var_grid[x]["Q"]
	    
	        # generating network
	        G = rmg.generate_modular_networks(n, sfunction, modfunction, q, m, d)
	        nx.write_graphml(G, "networks/Sah_try2_Q"+str(q)+"_N"+str(n)+"_rep"+str(y)+".graphml")

except Exception as ex: 
	pass 
