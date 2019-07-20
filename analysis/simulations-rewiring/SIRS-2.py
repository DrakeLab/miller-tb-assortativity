import networkx as nx
import random
import heapq
import scipy
import EoN
from collections import defaultdict
from collections import Counter
import matplotlib.pyplot as plt

N = 1000

# SPONTANEOUS transitions
H = nx.DiGraph()  
#I->R recover to R
H.add_edge('I.f', 'R.f', rate = 1.4)  # female  
H.add_edge('I.m', 'R.m', rate = 1.4)  # male

#R->S revert to S
H.add_edge('R.f', 'S.f', rate = 0.2)   # female
H.add_edge('R.m', 'S.m', rate = 0.2)   # male

# INDUCED transitions
J = nx.DiGraph()    
J.add_edge(('I.f', 'S.f'), ('I.f', 'I.f'), rate = .25)  # female infects female
J.add_edge(('I.m', 'S.m'), ('I.m', 'I.m'), rate = 1.75)  # male infects male  
J.add_edge(('I.m', 'S.f'), ('I.m', 'I.f'), rate = .25)  # male infects female   
J.add_edge(('I.f', 'S.m'), ('I.f', 'I.m'), rate = 1.75)  # female infects male   

G = nx.fast_gnp_random_graph(N, 5./(N-1))

# make a dictionary to start with 
IC = defaultdict(lambda: 'S.f')

for node in range(N/2, N):
    IC[node] = 'S.m'
    
for node in range(100):
    IC[node] = 'I.m'

for node in range(200,300):
    IC[node] = 'I.f'

return_statuses = ('S.f', 'S.m', 'I.f', 'I.m', 'R.f', 'R.m')

sim = EoN.Gillespie_Arbitrary(G, H, J, IC, return_statuses, tmax = 30)

###### Plot

t= sim[0]
S=sim[1] + sim[2]
I=sim[3] + sim[4]
R=sim[5] + sim[6]

plt.plot(t, I, label='Total Infecteds')
plt.plot(t, sim[3], label = 'Females')
plt.plot(t, sim[4], label = 'Males')

plt.legend()
plt.xlabel('$t$')
plt.ylabel('Infecteds')
plt.savefig('SIRS-test.png')

