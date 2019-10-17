###### SLIRS_tra.py simulates SLIRS model with differences btwn M&F transmission on SW & SF assorted networks ###### 
###### Current version takes roughly X hours ######

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

###### Model parameters ######

N = [1000]           # Network Size
R = [0, 0.3, 0.6, 0.9]         #, 0.6, 0.9 Assortativity coefficient (Newman)
Tau = [0.04, 0.075, 0.1]         #  S->L Baseline transmission rate 
Del = [100000, 1./10.]     # L->I Reactivation rate; 10000=>SIR, del~0=SLIR
Gam = 1./2.          # I->R Recovery rate
Psi = [0, 0.33]       # R->S Reversion rate; 0=SIR, sig>0=SIRS
i0 = 0.05            # proportion initially infected 
tsteps = 200         # set max time steps to run model for

# Male:female differences to explain male bias
Alph_t = [1.0, 1.5, 2.0]   # Ratio of male:female susceptibility

# Network types
nt = ["G", "SW"]

var_grid = list(ParameterGrid({'N' : N, 'R' : R, 'Tau': Tau,
                               'Psi' : Psi, 'Del' : Del,
                               'Alph_t': Alph_t, 'net_type' : nt}))

reps = 1 + 100 # Number of reps

for x in range(0, len(var_grid)):
    n=var_grid[x]["N"]
    r=var_grid[x]["R"]
    tau=var_grid[x]["Tau"]
    alph_t=var_grid[x]["Alph_t"]
    type_net = var_grid[x]["net_type"]

    delt=var_grid[x]["Del"]
    psi=var_grid[x]["Psi"]

    ###### Model transitions ######

    # SPONTANEOUS transitions H
    H = nx.DiGraph()  

    #L->I reactivation of latent infection
    H.add_edge('L.f', 'I.f', rate = delt)  # female  
    H.add_edge('L.m', 'I.m', rate = delt)  # male

    #I->R recover to R
    H.add_edge('I.f', 'R.f', rate = Gam)  # female  
    H.add_edge('I.m', 'R.m', rate = Gam)  # male

    #R->S revert to S
    H.add_edge('R.f', 'S.f', rate = psi)   # female
    H.add_edge('R.m', 'S.m', rate = psi)   # male

    #S->I spontaneous infection
    H.add_edge('S.f', 'I.f', rate = psi/50)   # female
    H.add_edge('S.m', 'I.m', rate = psi/50)   # male


    # INDUCED transitions
    J = nx.DiGraph()

    #S->L I infects S
    J.add_edge(('I.f', 'S.f'), ('I.f', 'L.f'), rate = (tau * 2.0) / (alph_t + 1.0))  # female infects female
    J.add_edge(('I.m', 'S.m'), ('I.m', 'L.m'), rate = (tau * 2.0 * alph_t) / (alph_t + 1.0))  # male infects male  
    J.add_edge(('I.m', 'S.f'), ('I.m', 'L.f'), rate = (tau * 2.0 * alph_t) / (alph_t + 1.0))  # male infects female   
    J.add_edge(('I.f', 'S.m'), ('I.f', 'L.m'), rate = (tau * 2.0) / (alph_t + 1.0))  # female infects male      
               

    for y in range(1, reps):
        
        ###### READ GRAPH ######
        
        G = nx.read_graphml(path="networks3/"+str(type_net)+"_"+str(r)+"N"+str(n)+"rep"+str(y)+".graphml")

        ###### SET INITIAL CONDITIONS ######
        # note: len(IC) needs to be = # of nodes
        
        IC = defaultdict(lambda: "S.f") # initialize all susceptible women

        for i in range(len(G)):
            if G.node["n"+str(i)]['sex'] == 1.0: # but if sex=1, they are male
                IC["n"+str(i)] = 'S.m'
                if np.random.uniform(0,1,1)[0] < i0:
                    IC["n"+str(i)] = 'I.m'
            else:
                IC["n"+str(i)] = 'S.f' 
                if np.random.uniform(0,1,1)[0] < i0:
                    IC["n"+str(i)] = 'I.f' # set some susceptible f to infected f

        # Set state variables to return
        return_statuses = ('S.f', 'S.m', 'L.f', 'L.m',
                           'I.f', 'I.m', 'R.f', 'R.m')

        ###### RUN SIMULATION ######
        
        sim = EoN.Gillespie_Arbitrary(G, H, J, IC, return_statuses, tmax = tsteps)

        ###### SAVE SIMULATION ######
        
        tots = sim
        with open("SLIRS/"+str(type_net)+"_R"+str(r)+"_tau"+str(tau)+"_del"+str(delt)+
                  "_alph_t"+str(alph_t)+
                  "_psi"+str(psi)+"_rep"+str(y)+".csv",'wb') as out:
            csv_out=csv.writer(out)
            csv_out.writerow(['t','S.f','S.m',
                              'L.f', 'L.m','I.f',
                              'I.m', 'R.f', 'R.m'])
            csv_out.writerows(zip(*tots))

