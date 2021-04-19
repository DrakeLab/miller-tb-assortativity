###### Parallelized version of SLIRS.py ######
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

def process_sah(f):
    n=f["N"]
    r=f["R"]
    tau=f["Tau"]
    type_net = f["net_type"]
    delt=f["Del"]
    psi=f["Psi"]
    omeg=f["Omeg"]
    y=f["rep"]
    alph=f["Alph_vals"]
    alph_type=f["Alph_types"]

    #return [n, r, tau, y]

    ###### Model transitions ######

    # SPONTANEOUS transitions H
    H = nx.DiGraph()

    #L->I progression of latent infection
    H.add_edge('L.f', 'I.f', rate = delt)  # female
    H.add_edge('L.m', 'I.m', rate = delt)  # male

    #R->I relapse of active infection 
    H.add_edge('R.f', 'I.f', rate = omeg)   # female
    H.add_edge('R.m', 'I.m', rate = omeg)   # male

    #R->S revert to S
    H.add_edge('R.f', 'S.f', rate = psi)   # female
    H.add_edge('R.m', 'S.m', rate = psi)   # male

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


    #G = nx.read_graphml(path="networks/sah_net"+str(r)+"_"+str(y)+".graphml")

    ###### make GRAPH ######
    failed="no"
    if r>0: 
        try:
            G = rmg.generate_modular_networks(n, sg.geometric_sequence, sg.regular_sequence, r, 2, 10)
        except (RuntimeError, TypeError, NameError):
            failed="yes"
            G = nx.gnp_random_graph(n, .1)
    else:
        G = nx.gnp_random_graph(n, .1)
        
    clus = nx.average_clustering(G)
    path_len = nx.average_shortest_path_length(G)
    #deg_mean = mean(nx.degree(G).values())
    deg_assort = nx.degree_assortativity_coefficient(G)

    ###### SET INITIAL CONDITIONS ######
    # note: len(IC) needs to be = # of nodes
    #G.nodes()

    IC = defaultdict(lambda: "S.m")
    st0=random.randrange(1,n-1) # 1 initially infected individual

    for i in range(len(G)):
        if i < 500:
            IC[i] = 'S.f' # first 500 nodes women
            if i==st0:
                IC[i] = 'I.f'
        else:
            IC[i] = 'S.m'
            if i==st0:
                IC[i] = 'I.m'

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

    m_rec_rat = end[8] # males rec
    f_rec_rat = end[7] # females rec
    m_inf_rat = end[6] # males inf
    f_inf_rat = end[5] # females inf
    lat = end[3] +end[4] # amount of latent at end of sim for SLIRS model

    results = [failed, n, r, tau, alph, alph_type, delt, psi, omeg, y,
               type_net, clus, path_len, deg_assort,
               peak, sim_dur, m_rec_rat,f_rec_rat, m_inf_rat,f_inf_rat, lat]
    return results

##### SET UP #####

## Model parameters ##
N = [1000]           # Network Size
R = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3] # Q=r/2 (Newman)
Tau = [0.01, 0.02, 0.03]         #  S->L Baseline transmission rate
Del = [100000, 1./10.]     # L->I Incubation; 10000=>SIR, del~0=SLIR
Gam = 1./6.          # I->R Recovery rate (6 months)
Omeg = [0, (1./24)*.27]  # 0 or 27% of recovered individuals relapse within 24 months
Psi =  [0, 1./(12*10)]       # R->S Reversion rate; 0=SIR, sig>0=SIRS
tsteps = 300         # set max time steps to run model for

# Male:female differences to explain male bias
Alph_vals =  [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5]   # Ratio of male:female susceptibility
Alph_types = ["SUS", "TRA", "INF_PER"]

# Network parameters
nt = ["SAH"]

reps = range(1,25) # Number of reps

var_grid = list(ParameterGrid({'N' : N, 'R' : R, 'Tau': Tau,
                               'Psi' : Psi, 'Omeg' : Omeg, 'Del' : Del,
                               'Alph_vals': Alph_vals,'Alph_types': Alph_types,
                               'net_type' : nt,
                               'rep': reps}))

p = multiprocessing.Pool(28) # create a pool of 2 workers

sim_results = p.map(process_sah, var_grid) # perform the calculations

#print(sim_results)

with open("SLIRS-res/sah_res_resubmit_test2.csv",'wb') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(["failed", "net_size", "r", "tau", "alph_val", "alph_type", "delt", "psi", "omeg", "rep","type_net", "net_clus", "net_path_len", "net_deg_assort", "peak", "sim_dur", "m_rec", "f_rec", "m_inf","f_inf", "tot_lat"])
    csv_out.writerows(sim_results)

def process_rewire(f):
    n=f["N"]
    r=f["R"]
    tau=f["Tau"]
    type_net = f["net_type"]
    delt=f["Del"]
    psi=f["Psi"]
    omeg=f["Omeg"]
    y=f["rep"]
    alph=f["Alph_vals"]
    alph_type=f["Alph_types"]


    ###### READ GRAPH ######
    failed="no"
    try:
        G = nx.read_graphml(path="rewired-networks/"+str(type_net)+"_"+str(r)+"N"+str(n)+"rep"+str(y)+".graphml")
    except (RuntimeError, TypeError, NameError, IOError):
        failed="yes"
        G = nx.gnp_random_graph(n, .1)
        results = [failed, n, r, tau, alph, alph_type, delt, psi, y,type_net, 0, 0, 0,0, 0, 0,0, 0,0, 0]
        return results

    clus = nx.average_clustering(G)
    path_len = nx.average_shortest_path_length(G)
    #deg_mean = mean(nx.degree(G).values())
    deg_assort = nx.degree_assortativity_coefficient(G)

    ###### Model transitions ######

    # SPONTANEOUS transitions H
    H = nx.DiGraph()

    #L->I progression of latent infection
    H.add_edge('L.f', 'I.f', rate = delt)  # female
    H.add_edge('L.m', 'I.m', rate = delt)  # male

    #R->S revert to S
    H.add_edge('R.f', 'S.f', rate = psi)   # female
    H.add_edge('R.m', 'S.m', rate = psi)   # male

    #R->I relapse of active infection
    H.add_edge('R.f', 'I.f', rate = omeg)   # female
    H.add_edge('R.m', 'I.m', rate = omeg)   # male

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
    st0=random.randrange(1,n-1) # 1 initially infected individual

    for i in range(len(G)):
        if G.node["n"+str(i)]['sex'] == 1.0: # if sex=1, they are male
            IC["n"+str(i)] = 'S.m'
            if i==st0:
                IC[i] = 'I.m'

        else:
            IC["n"+str(i)] = 'S.f'
            if i==st0:
                IC[i] = 'I.f'

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

    m_rec_rat = end[8] # males rec
    f_rec_rat = end[7] # females rec
    m_inf_rat = end[6] # males inf
    f_inf_rat = end[5] # females inf
    lat = end[3] +end[4] # amount of latent at end of sim for SLIRS model

    results = [failed, n, r, tau, alph, alph_type, delt, psi, omeg, y,
               type_net, clus, path_len, deg_assort,
               peak, sim_dur, m_rec_rat,f_rec_rat, m_inf_rat,f_inf_rat, lat]
    return results

##### SET UP REWIRE #####
# Model parameters ##
N = [1000]           # Network Size
R = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6] # Q=r/2 (Newman)
Tau = [0.01, .0.02, 0.03]         #  S->L Baseline transmission rate
Del = [100000, 1./10.]     # L->I Incubation; 10000=>SIR, del~0=SLIR
Gam = 1./6.          # I->R Recovery rate (6 months)
Omeg = [0, (1./24)*.27]  # 0 or 27% of recovered individuals relapse within 24 months
Psi =  [0, 1./(12*10)]       # R->S Reversion rate; 0=SIR, sig>0=SIRS
tsteps = 300         # set max time steps to run model for

# Male:female differences to explain male bias
Alph_vals =  [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5]   # Ratio of male:female susceptibility
Alph_types = ["SUS", "TRA", "INF_PER"]

# Network parameters
nt = ["G", "SW"]
reps = range(1,25)

var_grid = list(ParameterGrid({'N' : N, 'R' : R, 'Tau': Tau,
                               'Psi' : Psi, 'Omeg' : Omeg, 'Del' : Del,
                               'Alph_vals': Alph_vals,'Alph_types': Alph_types,
                               'net_type' : nt,
                               'rep': reps}))

p = multiprocessing.Pool(31) # create a pool of workers

sim_results = p.map(process_rewire, var_grid) # perform the calculations

#print(sim_results)

with open("SLIRS-res/rewired_res_resubmit_test2.csv",'wb') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(["failed", "net_size", "r", "tau", "alph_val", "alph_type", "delt", "psi", "omeg", "rep", "type_net", "net_clus", "net_path_len", "net_deg_assort", "peak", "sim_dur", "m_rec", "f_rec", "m_inf","f_inf", "tot_lat"])
    csv_out.writerows(sim_results)

