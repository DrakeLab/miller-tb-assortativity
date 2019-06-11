import networkx as nx
import EoN
import matplotlib.pyplot as plt
import numpy as np
import os
import csv

os.chdir('/Users/paigemiller/Documents/phd/research-projects/miller-tb-assortativity/analysis/modular_graph_generator-master')



G = nx.read_graphml(path="scalefree_Q0_N1000_d10_m2.graphml")

initial_size = 100
gamma = 1.
tau = 0.1
t, S, I, R = EoN.fast_SIR(G, tau,
                          gamma,
                          initial_infecteds = [str(x) for x in range(initial_size)])

plt.semilogy(t, S, label = 'Susceptible')
plt.semilogy(t, I, label = 'Infected')
plt.semilogy(t, R, label = 'Recovered')
plt.legend()
plt.savefig('SIR.png')
