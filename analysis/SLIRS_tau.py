###### SLIRS_tau.py simulates SLIRS model with varying baseline transmission rates to determine R0 across models ###### 
###### Current version takes roughly X hours ######

import random_modular_generator_variable_modules as rmg
import sequence_generator as sg
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
from statistics import mean 

def process_file(f):
    n=f["N"]
    r=f["R"]
    tau=f["Tau"]
    delt=f["Del"]
    psi=f["Psi"]
    y=f["rep"]

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
    H.add_edge('S.f', 'I.f', rate = psi/10)   # female
    H.add_edge('S.m', 'I.m', rate = psi/10)   # male

    # INDUCED transitions
    J = nx.DiGraph()

    #S->L I infects S
    J.add_edge(('I.f', 'S.f'), ('I.f', 'L.f'), rate = tau)  # female infects female
    J.add_edge(('I.m', 'S.m'), ('I.m', 'L.m'), rate = tau)  # male infects male  
    J.add_edge(('I.m', 'S.f'), ('I.m', 'L.f'), rate = tau)  # male infects female   
    J.add_edge(('I.f', 'S.m'), ('I.f', 'L.m'), rate = tau)  # female infects male      
    
    # ###### READ REWIRED GRAPH ######
    # G = nx.read_graphml(path="networks3/G_"+str(r)+"N"+str(n)+"rep"+str(y)+".graphml")
        
    ###### MAKE SAH GRAPHS ######
    failed="no"
    try: 
        G = rmg.generate_modular_networks(n, sg.geometric_sequence, sg.regular_sequence, r, 2, 10)
    except (RuntimeError, TypeError, NameError):
        failed="yes"
        G = nx.gnp_random_graph(n, .1) 
    
    ###### CALCULATE NETWORK STATS #####
    deg_dict = list(G.degree())
    deg_vals_1=[x[1] for x in deg_dict]
    net_mean_k=mean(deg_vals_1)   # <K>
    deg_vals_2=[x[1]**2 - x[1] for x in deg_dict]
    net_var_k=mean(deg_vals_2)    # <K^2-K>

    ###### SET INITIAL CONDITIONS ######
    # note: len(IC) needs to be = # of nodes
        
    IC = defaultdict(lambda: "S.m") 
    for i in range(len(G)):
        if i < 500: 
            IC[i] = 'S.f' # first 500 nodes women
            if np.random.uniform(0,1,1)[0] < i0:
                IC[i] = 'I.f'
        else:
            IC[i] = 'S.m' 
            if np.random.uniform(0,1,1)[0] < i0:
                IC[i] = 'I.m' # set some susceptible f to infected f

    # Set state variables to return
    return_statuses = ('S.f', 'S.m', 'L.f', 'L.m',
                           'I.f', 'I.m', 'R.f', 'R.m')

    ###### RUN SIMULATION ######
    sim = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = tsteps)

    ###### GET SIMULATION RESULTS ######
    end = zip(*sim)[-1] # ending values for SSLLIIRR
    sim_dur = end[0] # duration of simulation
    m_rec = end[8] # males rec
    f_rec = end[7] # females rec
    results = [failed, n, r, tau, delt, psi, y, net_mean_k, net_var_k, sim_dur, m_rec, f_rec]
    return results
        
###### Model parameters ######
N = [1000]           # Network Size
R = [0, 0.3]         # (modularity if using Sah)
Tau = [0.0001, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.06, 0.075, 0.1, 0.3]    #  S->L Baseline transmission rate 
Del = [100000, 1./4.]     # L->I Reactivation rate; 10000=>SIR, del~0=SLIR
Gam = 1./2.          # I->R Recovery rate
Psi = [0]       # R->S Reversion rate; 0=SIR, sig>0=SIRS
i0 = 0.01            # proportion initially infected 
tsteps = 200         # set max time steps to run model for

var_grid = list(ParameterGrid({'N' : N, 'R' : R, 'Tau': Tau,
                               'Psi' : Psi, 'Del' : Del}))

reps = range(0,26) # Number of reps

var_grid = list(ParameterGrid({'N' : N, 'R' : R, 'Tau': Tau,
                               'Psi' : Psi, 'Del' : Del,
                               'rep': reps}))

p = multiprocessing.Pool(31) # create a pool of 2 workers

sim_results = p.map(process_file, var_grid) # perform the calculations

with open("SLIRS-res/"+"sah_res_tau_resubmit_test.csv",'wb') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(['failed', 'n', 'q', 'tau', 'delt', 'psi', 'y','net_mean_k', 'net_var_k','sim_dur', 'm_rec', 'f_rec'])
    csv_out.writerows(sim_results)

