###### SIRS-1.py simulates SIRS model with variable susceptibility on assorted networks ###### 
###### Version 1 uses SEPARATE compartments for females and males ###### 
###### Version 2 uses Joel Miller's approach ###### 

import networkx as nx
import EoN
import random
import matplotlib.pyplot as plt
import heapq
import scipy
import numpy as np
import os
import csv
from sklearn.model_selection import ParameterGrid
from collections import defaultdict
from collections import Counter

os.chdir('/Users/paigemiller/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-rewiring')


###### Model parameters ######

N = [1000] # Network Size
R = [0, 0.3, 0.6, 0.9] # Assortativity coefficient (Newman)
Tau = [0.12, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36] #[0.08, 0.16, .24, .32] # Baseline transmission rate 
Gamma = 0.5 # Recovery rate
Alph = [1.0, 2.0, 3.0, 4.0, 5.0, 10.0] #[1.0, 1.5, 2.0, 2.5, 5.0, 10.0] # Ratio of male:female susceptibility
Sigma = 0.1 # Reversion to Susceptible
i0 = 0.05 # proportion initially infected 
tsteps = 200 # set max time steps to fun model for

var_grid = list(ParameterGrid({'N' : N, 'R' : R, 'Tau': Tau, 'Alph': Alph}))

reps = 300 +1 # Number of reps

for x in range(0, len(var_grid)):
    n=var_grid[x]["N"]
    r=var_grid[x]["R"]
    tau=var_grid[x]["Tau"]
    alph=var_grid[x]["Alph"]

    ###### Model transitions ######

    # SPONTANEOUS transitions H
    H = nx.DiGraph()  

    #I->R recover to R
    H.add_edge('I.f', 'R.f', rate = Gamma)  # female  
    H.add_edge('I.m', 'R.m', rate = Gamma)  # male

    #R->S revert to S
    H.add_edge('R.f', 'S.f', rate = Sigma)   # female
    H.add_edge('R.m', 'S.m', rate = Sigma)   # male

    # INDUCED transitionsJ
    J = nx.DiGraph()

    #S->I I infects S
    J.add_edge(('I.f', 'S.f'), ('I.f', 'I.f'), rate = (tau * 2.0) / (alph + 1.0))  # female infects female
    J.add_edge(('I.m', 'S.m'), ('I.m', 'I.m'), rate = (tau * 2.0 * alph) / (alph + 1.0))  # male infects male  
    J.add_edge(('I.m', 'S.f'), ('I.m', 'I.f'), rate = (tau * 2.0) / (alph + 1.0))  # male infects female   
    J.add_edge(('I.f', 'S.m'), ('I.f', 'I.m'), rate = (tau * 2.0 * alph) / (alph + 1.0))  # female infects male      
    
    for y in range(1, reps):
        
        ###### READ GRAPH ######

        G = nx.read_graphml(path="networks3/G_"+str(r)+"N"+str(n)+"rep"+str(y)+".graphml")

        ###### SET INITIAL CONDITIONS ######
        # note: len(IC) needs to be = # of nodes
        
        IC = defaultdict(lambda: "S.f") # all susceptible women

        for i in range(len(G)):
            if G.node["n"+str(i)]['sex'] == 1.0:
                IC["n"+str(i)] = 'S.m'
                if np.random.uniform(0,1,1)[0] < i0:
                    IC["n"+str(i)] = 'I.m'
            else:
                IC["n"+str(i)] = 'S.f' # some are susceptible men
                if np.random.uniform(0,1,1)[0] < i0:
                    IC["n"+str(i)] = 'I.f'

        # Set state variables to return
        return_statuses = ('S.f', 'S.m', 'I.f', 'I.m', 'R.f', 'R.m')

        ###### RUN SIMULATION ######
        
        sim = EoN.Gillespie_Arbitrary(G, H, J, IC, return_statuses, tmax = tsteps)

        ###### SAVE SIMULATION ######
        
        tots = sim
        with open("SIRS/SIRS1_R"+str(r)+"_N"+str(n)+"_tau"+str(tau)+"_alph"+str(alph)+"_rep"+str(y)+".csv",'wb') as out:
            csv_out=csv.writer(out)
            csv_out.writerow(['t','S.f','S.m', 'I.f','I.m', 'R.f', 'R.m'])
            csv_out.writerows(zip(*tots))



###### Test plot ######

##t= sim[0]
##S=sim[1] + sim[2]
##I=sim[3] + sim[4]
##R=sim[5] + sim[6]
##
##plt.plot(t, I, label='Total Infecteds')
##plt.plot(t, sim[3], label = 'Females')
##plt.plot(t, sim[4], label = 'Males')
##
##plt.legend()
##plt.xlabel('$t$')
##plt.ylabel('Infecteds')
##plt.savefig('SIRS-test2.png')


