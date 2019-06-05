### Import modules

import networkx as nx
import random_modular_generator_variable_modules as rmg
import sequence_generator as sg
import EoN
from sklearn.model_selection import ParameterGrid
import numpy as np
import os
import csv

os.chdir('/Users/paigemiller/Documents/phd/research-projects/miller-tb-assortativity/analysis/modular_graph_generator-master')

##### Parameters for network generation

# varying parameters
N= [5000] #[2000, 5000]                                                # network sizes
Q=  [-0.2, -0.1, 0, 0.1, 0.2, 0.4]            # modularity ~ assort/2  (already done)

var_grid= list(ParameterGrid({'N' : N, 'Q' : Q}))

# constant parameters (Network)
d= 10                                      # mean degree
m= 2                                       # number of modules (sexes)
sfunction = sg.scalefree_sequence          # degree distribution 
modfunction = sg.regular_sequence          # module size distribution
# loop variables
reps= 10

# constant parameters (SIR)
tau= 0.125           # transmission rate
gamma= 1.0         # recovery rate
rho= 0.05         # random fraction initially infected

for x in range(0, len(var_grid)): 
    for y in range(0, reps):
        n=var_grid[x]["N"]
        q=var_grid[x]["Q"]
    
        # generating network
        G = rmg.generate_modular_networks(n, sfunction, modfunction, q, m, d)
        nx.write_graphml(G, "G_Q"+str(q)+"_N"+str(n)+"_rep"+str(y)+".graphml")

        # running SIR model
        sim = EoN.fast_SIR(G, tau, gamma, rho=rho, return_full_data=True)
        tots = sim.summary()
        with open("SIR_Q"+str(q)+"_N"+str(n)+"_rep"+str(y)+".csv",'wb') as out:
            csv_out=csv.writer(out)
            csv_out.writerow(['t','s','i','r'])
            csv_out.writerows(zip(*tots))

        # extract data on node infection at last time step
        res = sim.get_statuses(time=sim.t()[-1])
        with open("Final_Q"+str(q)+"_N"+str(n)+"_rep"+str(y)+".csv",'wb') as csv_file:
           writer = csv.writer(csv_file)
           for key, value in res.items():
               writer.writerow([key, value])
