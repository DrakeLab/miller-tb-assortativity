###### Parallelized version of SLIRS.py ######

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
import multiprocessing 

def process_file(f):
    n=f["N"]
    r=f["R"]
    tau=f["Tau"]
    #type_net = f["net_type"]
    delt=f["Del"]
    psi=f["Psi"]
    y=f["rep"]
    alph=f["Alph_vals"]
    alph_type=f["Alph_types"]

    #return [n, r, tau, y]
    
    ###### READ GRAPH ######
    
    G = nx.read_graphml(path="~/Documents/miller-tb-assortativity/analysis/simulations-sah/networks/GG_Q"+str(r)+"_N1000"+"_rep"+str(y)+".graphml")

    clus = nx.average_clustering(G)
    path_len = nx.average_shortest_path_length(G)
    #deg_mean = 0 #mean(nx.degree(G).values())
    deg_assort = nx.degree_assortativity_coefficient(G)

    #return [clus]

    # ###### Model transitions ######

    # SPONTANEOUS transitions H
    H = nx.DiGraph()  

    #L->I reactivation of latent infection
    H.add_edge('L.f', 'I.f', rate = delt)  # female  
    H.add_edge('L.m', 'I.m', rate = delt)  # male

    #R->S revert to S
    H.add_edge('R.f', 'S.f', rate = psi)   # female
    H.add_edge('R.m', 'S.m', rate = psi)   # male

    #S->I spontaneous infection
    H.add_edge('S.f', 'I.f', rate = psi/50)   # female
    H.add_edge('S.m', 'I.m', rate = psi/50)   # male

    # INDUCED transitions
    J = nx.DiGraph()

    if alph_type == "SUS": 
        #S->L I infects S
        J.add_edge(('I.f', 'S.f'), ('I.f', 'L.f'), rate = (tau * 2.0) / (alph + 1.0))         # female infects female
        J.add_edge(('I.m', 'S.m'), ('I.m', 'L.m'), rate = (tau * 2.0 * alph) / (alph + 1.0))  # male infects male  
        J.add_edge(('I.m', 'S.f'), ('I.m', 'L.f'), rate = (tau * 2.0) / (alph + 1.0))         # male infects female   
        J.add_edge(('I.f', 'S.m'), ('I.f', 'L.m'), rate = (tau * 2.0 * alph) / (alph + 1.0))  # female infects male      

        #I->R recover to R
        H.add_edge('I.f', 'R.f', rate = Gam)  # female  
        H.add_edge('I.m', 'R.m', rate = Gam)  # male

    elif alph_type == "TRA":
        #S->L I infects S
        J.add_edge(('I.f', 'S.f'), ('I.f', 'L.f'), rate = (tau * 2.0) / (alph + 1.0))         # female infects female
        J.add_edge(('I.m', 'S.m'), ('I.m', 'L.m'), rate = (tau * 2.0 * alph) / (alph + 1.0))  # male infects male  
        J.add_edge(('I.m', 'S.f'), ('I.m', 'L.f'), rate = (tau * 2.0 * alph) / (alph + 1.0))  # male infects female   
        J.add_edge(('I.f', 'S.m'), ('I.f', 'L.m'), rate = (tau * 2.0) / (alph + 1.0))         # female infects male      

        #I->R recover to R
        H.add_edge('I.f', 'R.f', rate = Gam)  # female  
        H.add_edge('I.m', 'R.m', rate = Gam)  # male

    else: # alph_type=="INF_PER"
        #S->L I infects S
        J.add_edge(('I.f', 'S.f'), ('I.f', 'L.f'), rate = tau)  # female infects female
        J.add_edge(('I.m', 'S.m'), ('I.m', 'L.m'), rate = tau)  # male infects male  
        J.add_edge(('I.m', 'S.f'), ('I.m', 'L.f'), rate = tau)  # male infects female   
        J.add_edge(('I.f', 'S.m'), ('I.f', 'L.m'), rate = tau)  # female infects male      

        #I->R recover to R
        H.add_edge('I.f', 'R.f', rate = (Gam * (alph + 1))/2)         # female  
        H.add_edge('I.m', 'R.m', rate = (Gam * (alph + 1))/(2*alph))  # male
           

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
    
    sim = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = tsteps)

    ###### GET SIMULATION RESULTS ######
    tots = zip(*sim)

    s = []

    for row in tots:
        s.append(row[5] + row[6]) # time series of infecteds

    peak = max(s)
    
    end = zip(*sim)[-1] # ending values for SSLLIIRR

    sim_dur = end[0] # duration of simulation

    mf_rec_rat = end[8]/end[7] # MF ratio for SIR/SLIR models
    mf_inf_rat = end[6]/(end[5]+0.1) # MF ratio for SIRS/SLIRS models
    lat = end[3] +end[4] # amount of latent at end of sim for SLIRS model

    rec_size = end[7] + end[8] # for SIR/SLIR models
    inf_size = end[5] + end[6] # for SIRS/SLIRS models

    results = [n, r, tau, alph, alph_type, delt, psi, y,
               type_net, clus, path_len, deg_assort, 
               peak, sim_dur, mf_rec_rat, mf_inf_rat, lat, rec_size, inf_size]
    return results

# ##### SET UP #####

## Model parameters ##
N = [1000]           # Network Size
R = [0, 0.15, 0.3, 0.45]         #, 0.6, 0.9 Assortativity coefficient (Newman)
Tau = [0.04, 0.075, 0.1, 0.15]         #  S->L Baseline transmission rate 
Del = [100000, 1./10.]     # L->I Reactivation rate; 10000=>SIR, del~0=SLIR
Gam = 1./2.          # I->R Recovery rate
Psi =  [0, 0.33]       # R->S Reversion rate; 0=SIR, sig>0=SIRS
i0 = 0.05            # proportion initially infected 
tsteps = 250         # set max time steps to run model for

# Male:female differences to explain male bias
Alph_vals = [1.0,  2.0, 3.0]   # Ratio of male:female susceptibility
Alph_types = ["SUS", "TRA", "INF_PER"]

# Network parameters
#nt = ["G", "SW"]

reps = range(0,50) # Number of reps 

var_grid = list(ParameterGrid({'N' : N, 'R' : R, 'Tau': Tau,
                               'Psi' : Psi, 'Del' : Del,
                               'Alph_vals': Alph_vals,'Alph_types': Alph_types,
                               'rep': reps}))

p = multiprocessing.Pool(24) # create a pool of 2 workers

sim_results = p.map(process_file, var_grid) # perform the calculations

print(sim_results)

with open("test_res.csv",'wb') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(["n", "r", "tau", "alph_val", "alph_type", "reactivation_rate", "reversion_rate", "rep",
                      "type_net", "net_clustering", "net_path_len", "net_deg_assort", 
                      "peak", "duration", "mf_r_ratio", "mf_i_ratio", "latent_prev", "recovered_prev", "infected_prev"])
    csv_out.writerows(sim_results)
