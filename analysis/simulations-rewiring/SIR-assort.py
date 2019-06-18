import networkx as nx
import EoN
import matplotlib.pyplot as plt
import numpy as np
import os
import csv
from sklearn.model_selection import ParameterGrid

os.chdir('/Users/paigemiller/Documents/phd/research-projects/miller-tb-assortativity/analysis/network-generation/op-2-simulations/6-5-19')

###### Parameters for network generation

N = [500, 1000, 1500]
R = [-0.9, -0.8, -0.7, -0.6, -0.5,
     -0.4, -0.3, -0.2, -0.1,
     0,
      0.1,  0.2,  0.3,  0.4,  0.5,
      0.6,  0.7,  0.8,  0.9]
reps=100 +1

###### Parameters for epidemics on networks
Tau = [0.125, 0.25, 0.375]
Gamma = 1
Rho = 0.05

var_grid = list(ParameterGrid({'N' : N, 'R' : R, 'Tau': Tau}))

for x in range(0, len(var_grid)):
    for y in range(1, reps):

        n=var_grid[x]["N"]
        r=var_grid[x]["R"]
        tau=var_grid[x]["Tau"]

        G = nx.read_graphml(path="networks/G_"+str(r)+"N"+str(n)+"rep"+str(y)+".graphml")

        sim = EoN.fast_SIR(G, tau, Gamma, rho=Rho, return_full_data=True)

        tots = sim.summary()
        with open("epi-sims/SIR_R"+str(r)+"_N"+str(n)+"_tau"+str(tau)+"_rep"+str(y)+".csv",'wb') as out:
            csv_out=csv.writer(out)
            csv_out.writerow(['t','s','i','r'])
            csv_out.writerows(zip(*tots))

        # extract data on node infection at last time step
        res = sim.get_statuses(time=sim.t()[-1])
        with open("epi-sims/Final_R"+str(r)+"_N"+str(n)+"_tau"+str(tau)+"_rep"+str(y)+".csv",'wb') as csv_file:
           writer = csv.writer(csv_file)
           for key, value in res.items():
               writer.writerow([key, value])
