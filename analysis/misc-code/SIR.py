import networkx as nx
import EoN
import matplotlib.pyplot as plt

#G = nx.configuration_model([1,5,10]*100000)

#G = nx.barabasi_albert_graph(10000, 5)

G = nx.read_graphml(path="scalefree_Q0_N1000_d10_m2.graphml")

initial_size = 100
gamma = 1.
tau = 0.1
sigma=1/10



t, S, I, R = EoN.fast_SIR(G, tau,
                          gamma,
                          initial_infecteds = [str(x) for x in range(initial_size)])

plt.semilogy(t, S, label = 'Susceptible')
plt.semilogy(t, I, label = 'Infected')
plt.semilogy(t, R, label = 'Recovered')
plt.legend()
plt.savefig('SIR.png')
