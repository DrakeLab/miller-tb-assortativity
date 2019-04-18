---
header-includes:
   - \usepackage{physics}
output:
  pdf_document: default
  html_document: default
---

# PROTOCOL FOR: 

Sex-assortativity and the spread of TB in contact networks

$\textit{Spring 2019}$

## Authors: 

* Paige Miller
* John Drake
* TBD

### Background: 

* The global male:female TB case ratio is 1.9 and this pattern is strikingly consistent across countries
* This male-bias in TB case reports could be due to sex-specific differences in exposure to Mtb or susceptibility to disease following exposure, or a combination thereof
* Behavioral characteristics due to different societal gender roles may contribute to sex-specific rates of exposure to Mtb: differences in number (no evidence from Miller et al. _in prep_ or other studies Mossong et al. 2008, etc.) or types of contacts 
* Types of contacts could be mediated by assortative mixing which is where individuals tend to associate with others similar to them, creating sub-groups within a social network.
* In addition to exposure, behavioral and biological differences may contribute to sex-specific susceptibility to disease following exposure: smoking, immune responses, hormones
* We wondered whether the types of contacts that men have, specifically their preferential social mixing with other men, could contribute to male-bias in TB case reports per se or if previously described differences in susceptibility between the sexes are necessary to explain male-bias

### Research questions:

1. Can sex-assortative mixing lead to observed levels of overall TB prevalence and sex-specific TB prevalence?
2. Are sex-specific differences in susceptibility required to explain male-bias in TB cases? 

### Study design: 

We will examine the effects of sex-assortativity and sex-specific susceptibility on the ratio of male to female cases using models of disease spread on simulated social networks. 

| Variable  | Values tested  | 
|:-:|:-:|
| Assortativity, $r$  |  -0.4, - 0.2, 0, 0.2, 0.4 |
| Size, $N$ | 500, $1\cdot 10 ^ 3$  |

Table: Design of pilot study. 

_Network simulation option 1: de Almeida et al. (2013)_

We will generate scale-free networks that vary sex-assortativity using the algorithm of de Almeida et al. (2013). Their algorithm generates assorted networks with power-law degree distributions. The algorithm is as follows: 

1. Start with $m_0$ nodes which are characterized by some intrinsic characteristic $\eta \in (0, 1)$ which represents a node's hypothetical placement along a gender axis. Let $A_{ij} = | \eta_i - \eta_j |$ be the difference in intrinsic characteristics of nodes $i$ and $j$. 
2. At every time step, a new node $j$, with $\eta_j$ attributed randomly, attaches to other $m$ $(m \leq m_0)$ pre-existing nodes $i$ (note: unweighted, undirected, non-repeat edges) by considering jointly the connections between new sites and those having high number of nearest neighbors ($k_i$) and high similitude ($A_{ij}$) with probability 

\[ \pi_i = \frac{(1-A_{ij})k_i}{\sum_{i}(1-A_{ij})k_i} \]

3. This process is used to add the ($m_0 + 1$) site, ($m_0 + 2$), and so on until the network reaches the desired size, $N=m_0 + t$ where $t$ is the time variable.

This process of connecting nodes based on similarity in $\eta$ has not yet introduced "groups" or in our case, node sex. Following de Almeida et al. (2013), we assign nodes to be male (0) or female (1) nodes before the onset of adding nodes and before drawing edges based on $\eta$. Females are composed of all nodes have $\eta$ in the interval $\Delta \eta = \eta_{max} - \eta_{min}$ and males are composed of all nodes not in this interval (i.e. $\textbf {FILL IN EXPLANATION}$). 

Here, we will try $\textbf {INSERT VALUES FROM PAPER THAT MAKE SENSE}$ to get a feel for how values of these variables correspond to desired values for pilot study and then edit protocol accordingly. 

_Network simulation option 2: re-wire edges of basic SF networks_

We will generate scale-free networks that vary sex-assortativity using the classic BA-algorithm with the following adjustments: 

1. Assign nodes randomly as male (0) or female (1).  
2. Calculate network assortativity coefficient, r, from Newman (2003) where $\sum_{ij}{e_{ij}=1}$ is the sum of proportion of edges between nodes of each type (0 and 1), $\sum_{ij}{e_{ij}=a_i}$ is the proportion of 


__Simulate outbreaks on networks with and without differences in susceptibility by sex:__

* Gillespie like simulation on networks 
* Markov chain with exponential waiting times
* Doesn't need to have analytical solution

### Analysis: 

_Measure overall prevalence:_

_Measure ratio of male:female infections:_

### Checklist: 

* literature review on assortativity in social networks
* decide on how to model assorted networks
* decide on how to model epidemics on networks
* pilot study on small number or size of networks

### CHANGE-LOG:

* [bullet list of changes to protocols with reasons for each change]
