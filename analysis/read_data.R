# Data from epidemics on assorted networks
library(igraph)

########## TO DO ADD COLUMN TO CALCULATE R0 BEFORE NEXT RUN ######## see line 142 of analysis for example

###### ----------- SIR w/ constant susceptibility -----------  ###### 

ns <- c(5, 10, 15) * 100
rs <- 0:9 / 10
reps <- 1:100
tau <- c(.125, .25, .375)
vars <- expand.grid(r=rs, n=ns, rep=reps, tau=tau)

res <- data.frame(size=NA, rep=NA, tau=NA,
                  degAssort=NA, clustering=NA, pathLen=NA, 
                  r=NA, q=NA, 
                  time_steps=NA, peak_size=NA, 
                  peak_time=NA, tot_inf=NA, 
                  prev_1=NA, prev_2=NA, prev_ratio=NA)

sims <- list()

setwd("~/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-rewiring")

for(i in 1:nrow(vars)){
  
  s=vars[i, "n"]
  r=vars[i, "r"]
  rep=vars[i, "rep"]
  tau=vars[i, "tau"]
  
  ## Network data ###
  g <- read.graph(paste0("networks/G_",
                         r, "N",
                         s, "rep",
                         rep, ".graphml"),
                  format = "graphml")

  res[i, "size"] <- s
  res[i, "rep"] <- rep
  res[i, "tau"] <- tau
  res[i, "degAssort"] <- assortativity_degree(g, directed = FALSE)
  res[i, "clustering"] <- transitivity(g)
  res[i, "pathLen"] <- diameter(g, directed = FALSE)
  res[i, "r"] <- assortativity_nominal(g, types=V(g)$sex, directed=FALSE)
  res[i, "q"] <- modularity(g, membership=V(g)$sex)

  ### Epidemic data ### 
  # Combined state information (all time)
  simOut=read.csv(paste0("SIR/constant-suscept/SIR_R",vars[i, "r"], 
                         "_N", s, 
                         "_tau", tau, 
                         "_rep", vars[i, 'rep'], ".csv"))
  
  # simOut$tau <- tau
  # simOut$rep <- rep
  # simOut$r_lev <- r
  # sims[[i]] <- simOut
  
  # Final status of nodes 
  nodeOut=read.csv(paste0("SIR/constant-suscept/Final_R",vars[i, "r"],
                          "_N", s,
                          "_tau", tau,
                          "_rep", vars[i, 'rep'], ".csv"),
                   header = FALSE, col.names=c("node", "state"))

  nodeOut$node <- gsub("n", "", nodeOut$node)
  nodeOut$node <- as.numeric(nodeOut$node)
  
  nodeOut <- dplyr::arrange(nodeOut, node)
  nodeOut$sex <- V(g)$sex


  # Last value for state information (last time point for all ndoes)
  # simulation variables
  res[i, "time_steps"] <- tail(simOut$t, 1)
  res[i, "peak_size"] <- max(simOut$i)
  res[i, "peak_time"] <- simOut$t[which.max(simOut$i)]
  res[i, "tot_inf"] <- table(nodeOut[,2])[["R"]]/s

  res[i, "prev_1"] <- length(which(nodeOut$state[nodeOut$sex==1]=="R"))/(s/2)
  res[i, "prev_2"] <- length(which(nodeOut$state[nodeOut$sex==2]=="R"))/(s/2)
  res[i, "prev_ratio"] <- res[i, "prev_1"]/res[i, "prev_2"]
  
}

write.csv(res, "constant_results.csv")

###### ----------- SIR w/ variable susceptibility -----------  ###### 

# graphs from first set of simulations -- run if needed
rs <- 0:9/10
ns <- 1e3
reps <- 1:100
tau <- c(0.08, 0.12, 0.16, 0.2, .24, .28, .32, 0.36)
alpha <- c(1, 1.25, 1.5, 1.75)

vars <- expand.grid(r=rs, n=ns, rep=reps, tau=tau, alpha=alpha)

res <- data.frame(size=NA, rep=NA, tau=NA, alph=NA,
                  degAssort=NA, clustering=NA, pathLen=NA, 
                  r=NA, q=NA, 
                  time_steps=NA, peak_size=NA, 
                  peak_time=NA, tot_inf=NA, 
                  prev_1=NA, prev_2=NA, prev_ratio=NA)

sims <- list()

setwd("~/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-rewiring")

for(i in 1:nrow(vars)){
  
  s=vars[i, "n"]
  r=vars[i, "r"]
  rep=vars[i, "rep"]
  tau=vars[i, "tau"]
  alph=vars[i, "alpha"]
  
  ## Network data ###
  g <- read.graph(paste0("networks/G_",
                         r, "N",
                         s, "rep",
                         rep, ".graphml"),
                  format = "graphml")
  
  res[i, "r"] <- r
  res[i, "size"] <- s
  res[i, "rep"] <- rep
  res[i, "tau"] <- tau
  res[i, "alph"] <- alph
  res[i, "degAssort"] <- assortativity_degree(g, directed = FALSE)
  res[i, "clustering"] <- transitivity(g)
  res[i, "pathLen"] <- diameter(g, directed = FALSE)
  res[i, "r_actual"] <- assortativity_nominal(g, types=V(g)$sex, directed=FALSE)
  res[i, "q"] <- modularity(g, membership=V(g)$sex)
  
  ### Epidemic data ### 
  # Combined state information (all time)
  simOut=read.csv(paste0("SIR/var-suscept/SIR_R", r, 
                         "_N", s, 
                         "_tau", tau, "_alph", alph,
                         "_rep", rep, ".csv"))
  
  # Final status of nodes 
  nodeOut=read.csv(paste0("SIR/var-suscept/Final_R",vars[i, "r"],
                          "_N", s,
                          "_tau", tau, "_alph", alph,
                          "_rep", rep, ".csv"),
                   header = FALSE, col.names=c("node", "state"))
  
  nodeOut$node <- gsub("n", "", nodeOut$node)
  nodeOut$node <- as.numeric(nodeOut$node)
  
  nodeOut <- dplyr::arrange(nodeOut, node)
  nodeOut$sex <- V(g)$sex
  
  # Last value for state information (last time point for all ndoes)
  # simulation variables
  res[i, "time_steps"] <- tail(simOut$t, 1)
  res[i, "peak_size"] <- max(simOut$i)
  res[i, "peak_time"] <- simOut$t[which.max(simOut$i)]
  res[i, "tot_inf"] <- table(nodeOut[,2])[["R"]]/s
  
  res[i, "prev_1"] <- length(which(nodeOut$state[nodeOut$sex==1]=="R"))/(s/2)
  res[i, "prev_2"] <- length(which(nodeOut$state[nodeOut$sex==2]=="R"))/(s/2)
  res[i, "prev_ratio"] <- res[i, "prev_1"]/res[i, "prev_2"]
  
}

write.csv(res, "SIR/variable_results.csv") # has initial infected =1; other one has proportion of initial infecteds=0.05



###### ----------- SIR2 w/ variable susceptibility (MORE REPS) -----------  ###### 

# graphs from first set of simulations -- run if needed
rs <- c(0, 0.3, 0.6, 0.9)
ns <- 1e3
reps <- 1:300
tau <- c(0.08, 0.16, .24, .32)
alpha <- c(1, 1.25, 1.5, 1.75)

vars <- expand.grid(r=rs, n=ns, rep=reps, tau=tau, alpha=alpha)

res <- data.frame(size=NA, rep=NA, tau=NA, alph=NA,
                  degAssort=NA, clustering=NA, pathLen=NA, 
                  r=NA, q=NA, 
                  time_steps=NA, peak_size=NA, 
                  peak_time=NA, tot_inf=NA, 
                  prev_1=NA, prev_2=NA, prev_ratio=NA)

sims <- list()

setwd("~/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-rewiring")

for(i in 1:nrow(vars)){

  s=vars[i, "n"]
  r=vars[i, "r"]
  rep=vars[i, "rep"] 
  tau=vars[i, "tau"]
  alph=vars[i, "alpha"]
  
  repg <- rep %% 100 
  
  ## Network data ###
  g <- read.graph(paste0("networks/G_",
                         r, "N",
                         s, "rep",
                         repg, ".graphml"),
                  format = "graphml")
  
  res[i, "r"] <- r
  res[i, "size"] <- s
  res[i, "rep"] <- rep 
  res[i, "tau"] <- tau
  res[i, "alph"] <- alph
  res[i, "degAssort"] <- assortativity_degree(g, directed = FALSE)
  res[i, "clustering"] <- transitivity(g)
  res[i, "pathLen"] <- diameter(g, directed = FALSE)
  res[i, "r_actual"] <- assortativity_nominal(g, types=V(g)$sex, directed=FALSE)
  res[i, "q"] <- modularity(g, membership=V(g)$sex)
  
  
  ### Epidemic data ### 
  # Combined state information (all time)
  simOut=read.csv(paste0("SIR/var-suscept-2/SIR_R", r, 
                         "_N", s, 
                         "_tau", tau, "_alph", alph,
                         "_rep", rep, ".csv"))
  
  # Final status of nodes 
  nodeOut=read.csv(paste0("SIR/var-suscept-2/Final_R", vars[i, "r"],
                          "_N", s,
                          "_tau", tau, "_alph", alph,
                          "_rep", rep, ".csv"),
                   header = FALSE, col.names=c("node", "state"))
  
  nodeOut$node <- gsub("n", "", nodeOut$node)
  nodeOut$node <- as.numeric(nodeOut$node)
  
  nodeOut <- dplyr::arrange(nodeOut, node)
  nodeOut$sex <- V(g)$sex
  
  # Last value for state information (last time point for all ndoes)
  # simulation variables
  res[i, "time_steps"] <- tail(simOut$t, 1)
  res[i, "peak_size"] <- max(simOut$i)
  res[i, "peak_time"] <- simOut$t[which.max(simOut$i)]
  res[i, "tot_inf"] <- table(nodeOut[,2])[["R"]]/s
  
  res[i, "prev_1"] <- length(which(nodeOut$state[nodeOut$sex==1]=="R"))/(s/2)
  res[i, "prev_2"] <- length(which(nodeOut$state[nodeOut$sex==2]=="R"))/(s/2)
  res[i, "prev_ratio"] <- res[i, "prev_1"]/res[i, "prev_2"]
  
}

write.csv(res, "SIR/variable_results_2.csv") # has initial infected =1; reps =300





