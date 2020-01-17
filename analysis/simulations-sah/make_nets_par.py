### Import modules

import random_modular_generator_variable_modules as rmg
import sequence_generator as sg
import networkx as nx
import random
from sklearn.model_selection import ParameterGrid
import numpy as np
import heapq
import scipy
import os
import csv
from collections import defaultdict
from collections import Counter
import multiprocessing 

#os.chdir('/Users/paigemiller/Documents/UGA/phd/research-projects/miller-tb-assortativity/analysis/simulations-sah')


def process_file(f): # parallelized version of Sah et al. (PNAS 2019) network generation algorithm
	n=f["N"]
	q=f["Q"]
	y=f["Y"]

	G = rmg.generate_modular_networks(n, sg.geometric_sequence, sg.regular_sequence, q, 2, 10)
	nx.write_graphml(G, "networks/GG_Q"+str(q)+"_N"+str(n)+"_rep"+str(y)+".graphml")



##### Parameters for SIR simulations
# Network params
N= [1000]                                  # network sizes
Q= [0, 0.15, 0.3, 0.45]                    # network modularity 
Y= range(0, 150)

var_grid = list(ParameterGrid({'N' : N, 'Q' : Q, 'Y' : Y}))

##### Create and save networks  ########
p = multiprocessing.Pool(15) # create a pool of workers
sim_results = p.map(process_file, var_grid) # perform the calculations
