###### Parallelized SLIR simulator on CP networks ######
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


def process_rewire(f):
    n=f["N"]
    d=f["d"]
    k = f["k"]
    y=f["Rep"]
    tau=f["Tau"]
    sig=f["Sigma"]
    gam=f["Gamma"]


    ###### READ GRAPH ######
    failed="no"
    try:
        G = nx.read_graphml(path="nets/N_"+str(n)+"D_10d_"+str(d)+"k_"+str(k)+"rep_"+str(y)+".graphml")
    except (RuntimeError, TypeError, NameError, IOError):
        failed="yes"
        G = nx.gnp_random_graph(n, .1)
    deg_dict = list(G.degree())
    deg_vals_1=[x[1] for x in deg_dict]
    net_mean_k=mean(deg_vals_1)   # <K>
    deg_vals_2=[x[1]**2 - x[1] for x in deg_dict]
    net_var_k=mean(deg_vals_2)    # <K^2-K>

    clus = nx.average_clustering(G)
    path_len = nx.average_shortest_path_length(G)
    #deg_mean = mean(nx.degree(G).values())
    deg_assort = nx.degree_assortativity_coefficient(G)

    ###### Model transitions ######
    # SPONTANEOUS transitions H
    H = nx.DiGraph()

    #L->I reactivation of latent infection
    H.add_edge('L', 'I', rate = sig)  

    #I->R recover to R
    H.add_edge('I', 'R', rate = gam)  


    # INDUCED transitions
    J = nx.DiGraph()
    J.add_edge(('I', 'S'), ('I', 'L'), rate = tau)         

    ###### SET INITIAL CONDITIONS ######
    # note: len(IC) needs to be = # of nodes

    IC = defaultdict(lambda: 'S')
    for node in range(5):
        IC[node] = 'I'

    # Set state variables to return
    return_statuses = ('S', 'L', 'I', 'R')

    ###### RUN SIMULATION ######
    sim = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = tsteps)

    ###### GET SIMULATION RESULTS ######
    tots = zip(*sim)

    s = []
    for row in tots:
        s.append(row[4]) # time series of infecteds

    peak = max(s)
    end = zip(*sim)[-1] 
    sim_dur = end[0] # duration of simulation
    final_size = end[4] # amount of recovered at end of sim for SLIRS model

    results = [failed, n, d, k, y, net_var_k, net_mean_k, tau, sig, gam, peak, sim_dur, final_size]
    return results

##### SET UP REWIRE #####
## Model parameters ##
N = [50, 100, 500, 1000, 5000]           # Network Size
d = [0.25, 0.5, 0.75] 
K = [1,2,3,4]
Tau = [0.04, 0.075, 0.1]         #  S->L Baseline transmission rate 
Sigma = [1./6.]     # L->I Reactivation rate; 10000=>SIR, del~0=SLIR
Gamma = [1./6.]          # I->R Recovery rate
tsteps = 300         # set max time steps to run model for

# Network parameters
reps = range(1,3)

var_grid = list(ParameterGrid({'N' : N, 'd' : d, 
                               'K' : K, 'Tau': Tau, 'Sigma' : Sigma, 'Gamma' : Gamma,
                               'Rep': reps}))

p = multiprocessing.Pool(2) # create a pool of workers

sim_results = p.map(process_rewire, var_grid) # perform the calculations

#print(sim_results)

with open("cp_prelim_25.csv",'wb') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(["failed", "n", "d", "k", "y", "net_var_k", "net_mean_k", "tau", "sig", "gam", "peak", "sim_dur", "final_size"])
    csv_out.writerows(sim_results)
