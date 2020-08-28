# Data and Code to Reproduce Results in "Effects of assortative mixing and sex-traits on male-bias in tuberculosis: A modelling study" 

## Authors: 

* Paige B. Miller
* Chris C. Whalen
* John M. Drake

### Summary: 

Globally, Tuberculosis disease (TB) is more common among males than females. Recent research proposes
that differences in preferential social mixing by sex, or sex-assortativity, can alter infection patterns in TB. We
conducted a simulation study to see whether sex-assorted mixing patterns can explain the global ratio of
male:female TB cases and what factors might cause sex-disparities in infectious diseases to be sensitive to
assortative mixing. Simulations showed sex-assortativity alone cannot cause sex-bias in TB. However, we find
an effect of interaction between assortativity and sex-traits that suggests a role for behaviour to influence sex-
specific epidemiology of infectious diseases. In our study, the role of sex-assortativity was especially apparent
for slower spreading infectious diseases, like TB. We also examined how assortativity and sex-traits affect the
final outbreak size and other epidemic dynamics. These results are important for understanding when sex-
assortativity, a common feature across human populations, can change epidemiological patterns.

### Contents: 

* `main.pdf`: Current draft of manuscript
* `supp.pdf`: Current draft of supplementary material and figures
* `random_modular_generator_variable_modules.py`: code from Sah et al. 2014 [1] used to model assorted networks in main text
* `rewire_nets.R`: Rewiring algorithm used to generate assorted networks, results shown in supplementary materials
* `run_simulations.py`: main script used to simulate disease spread on Sah networks & rewired networks, extract simulation information, save results as csv's
* `SLIRS_tau.py`: script used to model outbreak size as a function of varying transmission rates (tau)
* `analysis/analysis.Rmd`: uses data from run_simulations.py and SLIRS_tau.py to generate figures in /male-bias-figs folder
* `analysis/rewired-networks`: simulated rewired networks used in simulations
* `analysis/sah-networks`: simulated sah networks used in simulations
* `analysis/SLIRS-res`: files used to generate figures in main and supplementary materials

### Instructions to reproduce results: 

* To rerun full set of network simulations in main text, run `rewire-nets.R` and `run_simulations.py`. Results were obtained with R Version 4.0.0 and python version 2 on a linux machine with 32 cores, 188G of memory, and a ATI Mobility Radeon HD 5430 graphics card. Users may wish to decrease the number of simulations and number of cores to produce a subset of results. 
* Figures can be reproduced without rerunning network simulations with data stored in `SLIRS-res` and code in `analysis.Rmd`. Required packages and versions used to produce results are listed at the top of `.Rmd` file. 

[1] Sah P, Singh LO, Clauset A, Bansal S. Exploring community structure in biological networks with random graphs. BMC Bioinformatics. BioMed Central; 2014;15(220).
