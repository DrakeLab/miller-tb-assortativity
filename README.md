---
header-includes:
- \usepackage{physics}
- \usepackage{sectsty} \sectionfont{\centering}
output: pdf_document
---

# PROTOCOL FOR: 
# Sex-assortative contact patterns and male-bias of Tuberculosis

_Last updated: August 26, 2019_

## Authors: 

* Paige Miller
* John Drake

### Background: 

TB is a respiratory-transmitted infectious disease that is heterogeneously distributed globally and concentrated in Africa and Asia. In addition, TB is heterogeneously distributed within populations; in particular, men face higher risk factors than women. The global male:female TB case ratio is 1.9, a pattern strikingly consistent across countries. 

Male-bias in TB case reports is thought to be due to sex-specific differences in exposure to Mtb or susceptibility to disease following exposure, or a combination of the two. Exposure rates could be mediated by differences in behaviors and societal gender roles. For example men may have more contacts or be more central in social networks but there is a lack of evidence for this (Mossong et al. 2008, etc.). Alternatively, sex-specific exposure rates could be driven by differences in the types of contacts that men compared to women. Specifically, exposure could be higher in men due to assortative mixing which is where individuals tend to associate with others similar to them, creating sub-groups within a social network. However, we lack studies analyzing the role of hypothesized exposure mechanisms at driving male-bias. 

The effects of assortativity (i.e., modularity) epidemic dynamics (e.g., total outbreak size or equilibrium prevalence) have been studied with math models. Results indicate that the effect of assortativity depends on assumptions about recovery and supply of susceptible individuals (SIR, SIS) as well as the distribution of subgroup size. In SIR systems, the total outbreak size of diseases with SIR dynamics decrease with increasing assortativity, perhaps due to a build-up of recovered nodes within-subgroups and transmission bottlenecks between subgroups (Nadini 2018). However, very high assortativity may be required to produce this "protective" effect of assortativity (Salathe & Jones 2015, Sah et al. 2017). In contrast, for SIS processes, infected nodes return to susceptible nodes and a build-up of recovered nodes does not occur, resulting in intensely connected subgroups which increase the equilibrium prevalence of infection (Nandini 2018). Subgroup size also affects disease spread: networks with larger subgroups increase epidemic size when compared with networks with constant subgroup size (Sah et al. 2017). The number of subgroups also affects how assortativity relates to outbreak size: outbreak size on less fragmented (fewer subgroups) networks are not affected by assortativity but outbreak size on more fragmented networks increases with assortativity (Sah et al. 2017). Most prior studies have focused on networks with multiple (in the range of 10 to 1000, Sah et al. 2017), highly connected subgroups. For our purposes, an open question remains about how disease spreads within lightly assorted networks with two subgroups (i.e., men and women) which potentially vary in their susceptibility to infection (i.e., are men "supersusceptibles", Kraft 2015). 

In this project, we will investigate whether sex-assortativity in social networks can contribute to male-bias or sex-specific differences in susceptibility to infection are necessary to explain male-bias. In addition, we will compare how assortativity changes epidemic trajectories in SIR and SIRS models. SIRS models have the ability to capture the persistance of TB in endemic areas. 

### Research questions:

1. Can sex-assortative mixing lead to male-bias by itself or is sex-specific susceptibility required to explain male-bias?
2. What are the effects of sex-assortative mixing on disease spread (peak size/time, variation in outbreak size/duration)? 

### Study design: 

We will examine the effects of sex-assortativity and sex-specific susceptibility on the ratio of male to female cases using SIR and SIRS models of disease spread on synthetic contact networks. In the synthetic contact networks, each node will represent an individual and each link is an interaction between individuals that infection can spread. 

_Network generation:_

To understand how increasing assortativity could change patterns of sex-bias in TB, we will will vary assortativity, r. Assortativity will be calculated according to Newman's discrete assortativity coefficient (Newman 2003) as  

\[ r = \frac{\sum_i{e_{ii}} - \sum_i{(a(i)^2)}}{1-\sum_i{(1-a(i)^2)}}  = \frac{Tr \textbf{E} - ||E^2||}{1 - ||E^2||}\]

where $e_{ij}$ is the proportion of edges connecting nodes in subgroup i to subgroup j (undirected), $a(i)=\sum_{j}{e_(ij)}$. Alternatively, if $\textbf{E}$ is the matrix of $e_{ij}$ and $||X||$ is the sum of all elements in a matrix $\textbf{X}$, then assortativity can be calculated with the proportion of edges within-sex ($Tr \textbf{E}$) and the proportion of edges that would be within groups if connections were random $||E^2||$. Here, if edges were distributed randomly, $Tr \textbf{E} = ||E^2|| = 0.5$ and $r=0$. Typically, r ranges from $-1 \gg r \geq 1$ because disassortative networks are much closer to random networks (r=0) than are assortative networks. 

Note: Network modularity (Q) is a similar measure of non-random mixing in networks: $Q=\sum_{i}{e_{ii}-a^2_i}=Tr \textbf{E} - ||E^2||$ where $e_{ij}$ is the proportion of edges in the network that link nodes in community i to community j and $a_i=\sum_{j}{e_{ij}}$ represents the proportion of edges in the network that link to nodes in subgroup i. The maximum value of Q is $1-\frac{1}{K}$ where K is the number of modules (Sah 2014). Thus, the relationship between assortativity and modularity for networks with two subgroups is 

\[ assortativity = modularity / (1 - expected.prop.within.edges) \]

where the expected proportion of within edges is the proportion of nodes in that subgroup. Since groups have equal size, you would expect 1/2 to occur within group by chance. Thus, assortativity gets divided by 0.5 while modularity does not. 

We will generate scale-free networks according to the parameters listed in Table 1 using the classic BA-algorithm. Following network generation, we will update the networks as following: 

1. Assign nodes randomly as male (0) or female (1). 
2. Calculate temporary value of sex-assortativity in the network ($r_t$). 
3. If $r_t$ is not within $\epsilon$ of $r$, randomly choose a proportion $\alpha$ of 0--1 edges (i.e, a male--female edge) if $r\geq 0$ and re-wire them or if $r < 0$, choose a proportion $\alpha$ of 0--0 and 1--1 edges and re-wire them. 
4. Repeat step 3 until $|r_f - r_t| \leq \epsilon$. 

Network will be generated with parameters shown in Table 1.

| Variable  | Value  | 
|:-:|:-:|
| Sex-assortativity, $r$  | (0, 0.9) by 0.3 |
| Degree distribution, $p(K)$ | $\frac{k^{-\alpha}}{\zeta (\alpha)}$  |
| Mean degree, $<K>$ | 10  |
| Network size, $N$ | $1\cdot 10 ^ 3$|
| Tolerance, $\epsilon$ | 0.035  |
| Rewiring proportion, $\alpha$ | 0.2  |
| Replicates | 300  |

Table: Network parameters. Parameter range will be extended for publication.    

_Models of disease transmission:_ 

An infection spreads along an edge at probability depending on the baseline transmission rate, $\tau$, the susceptibility of the target node, and the duration of infection, T: $p(T)=1-e^{- \tau \sigma T}$. Infected individuals recover at an exponentially distributed recovery rate $\gamma$ and revert to susceptible at exponentially distributed reversion rate $\delta$ (set to 0 for SIR model). Parameter values for models are given in Table 2. 

| Variable  | Value  | 
|:-:|:-:|
| Initial susceptible, $S_0$  |  $N - 0.05 \cdot N$ |
| Initial infected, $I_0$  |  $0.05 \cdot N$ |
| Infection rate, $\tau$  |  0.025, 0.045, 0.08, 0.12 |
| Recovery rate, $\gamma$  |  1 |
| Reversion rate, $\delta$  |  0, 0.1 |
| M:F susceptibility ratio, $\alpha$  |  1, 2, 3 |

Table: Disease parameters for models. Parameter range will be extended for publication. 

Models will assume a non-Markovian process at first because the transmission rate will vary depending on susceptibility of each node. Once susceptibility is assigned, the process can be treated as Markovian (see Kiss, Miller, and Simon, page 224) and will be implemented with an event-based algorithm. Node susceptibility $\sigma$ will alter the transmission rate towards the target node as follows: 

\[ \sigma_f = \frac{2}{\alpha + 1} \]
\[ \sigma_m = \frac{2 \alpha}{\alpha + 1} \]

where f and m represent female and male nodes and $\alpha$ represents the susceptibility ratio between males and females. These equations satisfy 

\[ 0.5 \sigma_m + 0.5 \sigma_f = 1 \]
\[ \sigma_m = \alpha \sigma_f \]

so that the average susceptibility modifier is 1 and is independent of $R_0$ (Kiss, Miller, & Simon 2017). 

Parameters for $\tau$ were identified with the following analytical expression for $R_0$ in a continuous-time Markov model for SIR diseases (see Kiss, Miller, and Simon, page 221): 

\[ \frac{\tau}{\tau + \gamma} \frac{<K^2 - K>}{<K>} > 1  \]

Numerical validation for the solution will be obtained by simulating model with transmission rates leading above and below $R_0=1$. 

### Analysis: 

_Network structure:_

* Validation of network structure will include checking mean degree, number of edges, number of nodes, connectivity, presence of multiple edges, self-loops, etc. 
* Relationship between network structure (clustering, degree assortativity, diameter, degree variation) and assortativity

_M:F ratio:_

* Interaction between assortativity and susceptibility on M:F case ratio
* Variation by model type (SIR, SIRS)
* Variation in transmission rate

_Disease spread:_

* Relationship between assortativity and outbreak size, duration (for SIR model)
* Variation by network size and transmission rate

### Checklist: 

* X Understand the relationship between measures of community structure (Q vs. r) for K=2 modules
* X Run study across larger parameter grid and more replicates
* X Analyze and interpret results from extended simulations
* X Decide on next steps which could be: (1) Seed epidemics disproportionately in one module; (2) Incorporate sex-specific susceptibility; (3) Incorporate latent class of individuals; (4) Sample epidemics according to COHSONET and validate results
* X Edit protocol to add variable susceptibility for SIR
* X Write script for SIR model with variable susceptibility 
* X Run extended analysis of variable susceptibility for SIR (n=100)
* X Write script for SIRS model with variable susceptibility 
* X Run pilot of SIRS model with variable susceptibility (n=100)
* X Run analysis of SIR model with limited parameter range and more reps (n=300)
* X Write script for SIR and SIRS models for more comparative analysis (n=300)
* X Analyze results for SIR and SIRS comparative analysis
* X Write up analysis with variable susceptibility
* Understand and relate results to Salathe and Sah research 

### Important background papers: 

Kiss, I Z, J C Miller, and PL Simon Cham Springer. 2017. “Mathematics of Epidemics on Networks.” Springer.

Miller, Joel C. 2007. “Epidemic Size and Probability in Populations with Heterogeneous Infectivity and Susceptibility..” Physical Review. E, Statistical, Nonlinear, and Soft Matter Physics 76 (1 Pt 1): 010101. doi:10.1103/PhysRevE.76.010101.

Newman, MEJ. 2003. “Mixing Patterns in Networks.” Physical Review E 67 (2). American Physical Society. doi:10.1103/PhysRevE.67.026126.

Nhamoyebonde, Shepherd, and Alasdair Leslie. 2014. “Biological Differences Between the Sexes and Susceptibility to Tuberculosis.” Journal of Infectious Diseases 209 (suppl 3): S100–S106. doi:10.1093/infdis/jiu147.

Pastor-Satorras, R, and A Vespignani. 2001. “Epidemic Dynamics and Endemic States in Complex Networks.” Physical Review. E, Statistical, Nonlinear, and Soft Matter Physics 63 (6 Pt 2): 066117. doi:10.1103/PhysRevE.63.066117.

Perkins, S E, M F Ferrari, and P J Hudson. 2008. “The Effects of Social Structure and Sex-Biased Transmission on Macroparasite Infection.” Parasitology 135 (13): 1561–69. doi:10.1017/S0031182008000449.

Sah, Pratha, Stephan T Leu, Paul C Cross, Peter J Hudson, and Shweta Bansal. 2017. “Unraveling the Disease Consequences and Mechanisms of Modular Structure in Animal Social Networks..” Proceedings of the National Academy of Sciences of the United States of America 114 (16). National Academy of Sciences: 4165–70. doi:10.1073/pnas.1613616114.

Sah, Pratha, Lisa O Singh, Aaron Clauset, and Shweta Bansal. 2014. “Exploring Community Structure in Biological Networks with Random Graphs.” BMC Bioinformatics 15 (220). BioMed Central. doi:10.1186/1471-2105-15-220.

Salathé, Marcel, and James H Jones. 2010. “Dynamics and Control of Diseases in Networks with Community Structure.” Edited by Christophe Fraser. PLoS Computational Biology 6 (4). Public Library of Science: e1000736. doi:10.1371/journal.pcbi.1000736.

### CHANGE-LOG:

* Found little variation in M:F case ratio with assortativity only so next stage will compare the effects of sex-specific susceptibility and assortativity 
* SIR model on networks died out fairly quickly, not representative of TB in populations so next stage will compare SIR and SIRS models on networks (July 6)
* Found little variation in disease spread on networks of different sizes so next stage will focus on networks with 1000 nodes (July 6)
* Found results to be sensitive to transmission rates, next stages will focus on epidemics with varying transmission rates (July 24)
