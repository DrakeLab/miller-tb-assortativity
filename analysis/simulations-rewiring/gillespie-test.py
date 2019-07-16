import networkx as nx
import random
import heapq
import scipy
import EoN
from collections import defaultdict
from collections import Counter
import matplotlib.pyplot as plt

class _ListDict_(object):
    r'''
    The Gillespie algorithm with rejection-sampling will involve a step 
    that samples a random element from a set.  This is slow in Python.  
    So I'm introducing a new class based on a stack overflow answer by
    Amber (http://stackoverflow.com/users/148870/amber) 
    for a question by
    tba (http://stackoverflow.com/users/46521/tba) 
    found at
    http://stackoverflow.com/a/15993515/2966723
    '''
    def __init__(self, weighted = False):
        self.item_to_position = {}
        self.items = []

        self.weighted = weighted
        if self.weighted:
            self.weight = defaultdict(int) #presume all weights positive
            self.max_weight = 0
            self._total_weight = 0
            self.max_weight_count = 0
            

    def __len__(self):
        return len(self.items)

    def __contains__(self, item):
        return item in self.item_to_position

    def _update_max_weight(self):
        C = Counter(self.weight.values())  #may be a faster way to do this, we only need to count the max.
        self.max_weight = max(C.keys())
        self.max_weight_count = C[self.max_weight]

        
    def insert(self, item, weight = None):
        r'''
        If not present, then inserts the thing (with weight if appropriate)
        if already there, replaces the weight unless weight is 0
        
        If weight is 0, then it removes the item and doesn't replace.
        
        WARNING:
            replaces weight if already present, does not increment weight.
            
        
        '''
        if self.__contains__(item):
            self.remove(item)
        if weight != 0:
            self.update(item, weight_increment=weight)
        

    def update(self, item, weight_increment = None):
        r'''
        If not present, then inserts the thing (with weight if appropriate)
        if already there, increments weight
        
        WARNING:
            increments weight if already present, cannot overwrite weight.
        '''
        if weight_increment is not None: #will break if passing a weight to unweighted case
            if weight_increment >0 or self.weight[item] != self.max_weight:
                self.weight[item] = self.weight[item] + weight_increment
                self._total_weight += weight_increment
                if self.weight[item] > self.max_weight:
                    self.max_weight_count = 1
                    self.max_weight = self.weight[item]
                elif self.weight[item] == self.max_weight:
                    self.max_weight_count += 1
            else: #it's a negative increment and was at max
                self.max_weight_count -= 1
                self.weight[item] = self.weight[item] + weight_increment
                self._total_weight += weight_increment
                self.max_weight_count -= 1 
                if self.max_weight_count == 0:
                    self._update_max_weight               
        elif self.weighted:
            raise Exception('if weighted, must assign weight_increment')

        if item in self: #we've already got it, do nothing else
            return
        self.items.append(item)
        self.item_to_position[item] = len(self.items)-1

    def remove(self, choice):
        position = self.item_to_position.pop(choice)
        last_item = self.items.pop()
        if position != len(self.items):
            self.items[position] = last_item
            self.item_to_position[last_item] = position
            
        if self.weighted:
            weight = self.weight.pop(choice)
            self._total_weight -= weight
            if weight == self.max_weight:  
                #if we find ourselves in this case often
                #it may be better just to let max_weight be the
                #largest weight *ever* encountered, even if all remaining weights are less
                #
                self.max_weight_count -= 1
                if self.max_weight_count == 0 and len(self)>0:
                    self._update_max_weight()

    def choose_random(self):
        # r'''chooses a random node.  If there is a weight, it will use rejection
        # sampling to choose a random node until it succeeds'''
        if self.weighted:
            while True:
                choice = random.choice(self.items)
                if random.random() < self.weight[choice]/self.max_weight:
                    break
            # r = random.random()*self.total_weight
            # for item in self.items:
            #     r-= self.weight[item]
            #     if r<0:
            #         break
            return choice

        else:
            return random.choice(self.items)
        

    def random_removal(self):
        r'''uses other class methods to choose and then remove a random node'''
        choice = self.choose_random()
        self.remove(choice)
        return choice

    def total_weight(self):
        if self.weighted:
            return self._total_weight
        else:
            return len(self)
    def update_total_weight(self):
        self._total_weight = sum(self.weight[item] for item in self.items)



N = 100000
G = nx.fast_gnp_random_graph(N, 5./(N-1))

for node in G:
    if node < (N/2):
        G.node[node]["sex"] = "male"
        G.node[node]["sus"] = (2 * 1.5) / (1.5 + 1)

    else:
        G.node[node]["sex"] = "female"
        G.node[node]["sus"] = 2 / (1.5 + 1)

node_attribute_dict = {node: G.node[node]['sus'] for node in G.nodes()}

nx.set_node_attributes(G, values=node_attribute_dict, name='suscept_weight')

H = nx.DiGraph()  #DiGraph showing possible transitions that don't require an interaction
H.add_edge('I', 'R', rate = 1.4)   #I->R
H.add_edge('R', 'S', rate = 0.2)   #R->S

J = nx.DiGraph()    #DiGraph showing transition that does require an interaction.
J.add_edge(('I', 'S'), ('I', 'I'), rate = 1)  #IS->II

IC = defaultdict(lambda: 'S')
for node in range(200):
    IC[node] = 'I'
        
return_statuses = ('S', 'I', 'R')

tmin=0
tmax=10
spontaneous_transition_graph=H
nbr_induced_transition_graph=J

status = {node: IC[node] for node in G.nodes()}

    
times = [tmin]
data = {}
C = Counter(status.values())

for return_status in return_statuses:
        data[return_status] = [C[return_status]] 

spontaneous_transitions = list(spontaneous_transition_graph.edges())

induced_transitions = list(nbr_induced_transition_graph.edges())
potential_transitions = {}
rate = {}# intrinsic rate of a transition

# XXXXX dont understand get_weight function 
get_weight = defaultdict(lambda: defaultdict(lambda:None))


# i'm not giving weight labels to spontaneous transitions so got rid of it
for transition in spontaneous_transitions:
    rate[transition] = spontaneous_transition_graph.adj[transition[0]][transition[1]]['rate']
    potential_transitions[transition] = _ListDict_() #max_weight[transition]=1

# defining induced_transitions FIRST NODE KEEPS SAME STATUS
for transition in induced_transitions:
    if transition[0][0] != transition[1][0]:
        raise EoN.EoNError("transition {} -> {} not allowed: first node must keep same status".format(transition[0],transition[1]))
    rate[transition] = nbr_induced_transition_graph.adj[transition[0]][transition[1]]['rate'] 

    potential_transitions[transition] = _ListDict_()

#initialize all potential events to start.                
for node in G.nodes():        
    if spontaneous_transition_graph.has_node(status[node]) and spontaneous_transition_graph.degree(status[node])>0:
        for transition in spontaneous_transition_graph.edges(status[node]):
            potential_transitions[transition].update(node, weight_increment = get_weight[transition][node])
            #weight increment defaults to None if not present                    
                
                
    for nbr in G.neighbors(node):

        #print(status[node],status[nbr])
        if nbr_induced_transition_graph.has_node((status[node],status[nbr])) and nbr_induced_transition_graph.degree((status[node],status[nbr])) >0:
            for transition in nbr_induced_transition_graph.edges((status[node],status[nbr])):
                if (node, nbr) not in get_weight[transition]: #since edge may be in opposite order to earlier
                    get_weight[transition][(node, nbr)] = get_weight[transition][(nbr, node)]
                potential_transitions[transition].update((node, nbr), weight_increment = get_weight[transition][(node, nbr)])

t = tmin
