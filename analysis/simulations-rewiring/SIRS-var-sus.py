import networkx as nx
import EoN
import matplotlib.pyplot as plt
import numpy as np
import os
import csv
from sklearn.model_selection import ParameterGrid
from collections import defaultdict


os.chdir('/Users/paigemiller/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-rewiring')

###### Parameters for network generation

N = [50000] #[500, 1000, 1500]
##R = [-0.9, -0.8, -0.7, -0.6, -0.5,
##     -0.4, -0.3, -0.2, -0.1,
##     0,
##      0.1,  0.2,  0.3,  0.4,  0.5,
##      0.6,  0.7,  0.8,  0.9]
R = [0]  # Network assortativity (Can be one of above)

###### Parameters for epidemics on networks

Rho = [0.05]  # Determine initial infected proportion
Tau = [0.35] # IS->II
Gamma = [1.0]   # I->R
Sigma = [0.2] # R->S

var_grid = list(ParameterGrid({'N' : N, 'R' : R, 'Tau': Tau, 'Sigma': Sigma,
                               'Gamma': Gamma, 'Rho': Rho}))

reps=1 #100 +1

for x in range(0, len(var_grid)):
    for y in range(1,reps+1):

        n=var_grid[x]["N"]
        r=var_grid[x]["R"]
        tau=var_grid[x]["Tau"]
        gamma=var_grid[x]["Gamma"]
        sigma=var_grid[x]["Sigma"]
        i0=var_grid[x]["Rho"] * n

        #G = nx.read_graphml(path="networks/G_"+str(r)+"N"+str(n)+"rep"+str(y)+".graphml")
        G = nx.fast_gnp_random_graph(n, 5./(n-1))

        H = nx.DiGraph()  #DiGraph showing possible transitions that don't require an interaction
        H.add_edge('I', 'R', rate = gamma)   #I->R
        H.add_edge('R', 'S', rate = sigma)   #R->S

        J = nx.DiGraph()    #DiGraph showing transition that does require an interaction.
        J.add_edge(('I', 'S'), ('I', 'I'), rate = tau)  #IS->II

        IC = defaultdict(lambda: 'S')
        for node in range(int(i0)):
            IC[node] = 'I'

        return_statuses = ('S', 'I', 'R')

        t, S, I, R = EoN.Gillespie_Arbitrary(G, H, J, IC, return_statuses, tmax = 30)

        plt.plot(t, S, label = 'Susceptible')
        plt.plot(t, I, label = 'Infected')
        plt.plot(t, R, label = 'Recovered')
        plt.legend()
        plt.savefig('SIRS/SIRS' + str(y) + '.png')
        plt.clf()

##        tots = sim.summary()
##        with open("SIRS/SIRS_R"+str(r)+"_N"+str(n)+"_tau"+str(tau)+"_rep"+str(y)+".csv",'wb') as out:
##            csv_out=csv.writer(out)
##            csv_out.writerow(['t','s','i','r'])
##            csv_out.writerows(zip(*tots))
##
##        # extract data on node infection at last time step
##        res = sim.get_statuses(time=sim.t()[-1])
##        with open("SIRS/Final_R"+str(r)+"_N"+str(n)+"_tau"+str(tau)+"_rep"+str(y)+".csv",'wb') as csv_file:
##           writer = csv.writer(csv_file)
##           for key, value in res.items():
##               writer.writerow([key, value])
