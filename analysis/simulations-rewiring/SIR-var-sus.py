###### SIR-var-sus.py simulates SIR model with variable susceptibility on assorted networks ###### 

import networkx as nx
import EoN
import random
import matplotlib.pyplot as plt
import numpy as np
import os
import csv
from sklearn.model_selection import ParameterGrid

os.chdir('/Users/paigemiller/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-rewiring')

###### Model parameters ######

N = [1000] # Network Size
R = [0, 0.1, 0.2, 0.3,  0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

Tau = [0.08, 0.12] #[0.16, 0.2, .24, .28, .32, 0.36] # Baseline transmission rate 
Gamma = 1 # Recovery rate
Alph = [1, 1.25, 1.5, 1.75, 2] # Ratio of male:female susceptibility

var_grid = list(ParameterGrid({'N' : N, 'R' : R, 'Tau': Tau, 'Alph': Alph}))

reps = 100 + 1 # Number of reps

##### Model functions ########

def transmission_rate(source, target, tau):
    rate =  tau * G.node[target]['sus']
    return random.expovariate(rate)

def rec_time_function(node, gamma):
    return random.expovariate(gamma)

############# RUN MODEL, SAVE RESULTS #################

for x in range(0, len(var_grid)):
    for y in range(1, reps):

        n=var_grid[x]["N"]
        r=var_grid[x]["R"]
        tau=var_grid[x]["Tau"]
        alph=var_grid[x]["Alph"]

        G = nx.read_graphml(path="networks/G_"+str(r)+"N"+str(n)+"rep"+str(y)+".graphml")

        for node in G:
            if G.node[node]['sex'] == 1.0:
                G.node[node]['sus'] = (2 * alph) / (alph + 1)

            else:
                G.node[node]['sus'] = 2 / (alph + 1)


        sim = EoN.fast_nonMarkov_SIR(G, transmission_rate, rec_time_function,
                                     trans_time_args=(tau,), rec_time_args=(Gamma,),
                                     return_full_data=True)

        tots = sim.summary()
        with open("SIR/var-suscept/SIR_R"+str(r)+"_N"+str(n)+"_tau"+str(tau)+"_alph"+str(alph)+"_rep"+str(y)+".csv",'wb') as out:
            csv_out=csv.writer(out)
            csv_out.writerow(['t','s','i','r'])
            csv_out.writerows(zip(*tots))

        # extract data on node infection at last time step
        res = sim.get_statuses(time=sim.t()[-1])
        with open("SIR/var-suscept/Final_R"+str(r)+"_N"+str(n)+"_tau"+str(tau)+"_alph"+str(alph)+"_rep"+str(y)+".csv",'wb') as csv_file:
           writer = csv.writer(csv_file)
           for key, value in res.items():
               writer.writerow([key, value])

##### WORKING EXAMPLE OF SIR WITH SUSCEPTIBILITY ######


##os.chdir('/Users/paigemiller/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-rewiring')
##
##r=0.9
##n=1500
##y=21
##
##tau=0.4
##alph=2
##gamma=1
##
##G = nx.read_graphml(path="networks/G_"+str(r)+"N"+str(n)+"rep"+str(y)+".graphml")
##
##for node in G:
##    if G.node[node]['sex'] == 1.0:
##        G.node[node]['sus'] = tau * alph
##    else:
##        G.node[node]['sus'] = tau 
##
##def transmission_rate(source, target, tau):
##    rate =  tau * G.node[target]['sus']
##    return random.expovariate(rate)
##
##def rec_time_function(node, gamma):
##    return random.expovariate(gamma)
##
##sim = EoN.fast_nonMarkov_SIR(G, transmission_rate, rec_time_function,
##                            trans_time_args=(tau,), rec_time_args=(gamma,),
##                            rho = 0.01, return_full_data=True)
##
##
##
##t, S, I, R = sim.summary()
##plt.plot(t, I, label='Total Infecteds')
##
##t1, S1, I1, R1 = sim.summary(nodelist = [node for node in G if G.node[node]['sex']==1.0])
##plt.plot(t1, I1, label = 'MEN')
##
##t2, S2, I2, R2 = sim.summary(nodelist = [node for node in G if G.node[node]['sex']==2.0])
##plt.plot(t2, I2, label = 'WOMEN')
##
##plt.legend()
##plt.xlabel('$t$')
##plt.ylabel('Infecteds')
##plt.savefig('bipartite.png')
##
