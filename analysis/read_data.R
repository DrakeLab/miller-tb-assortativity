# Data from epidemics on assorted networks
library(igraph)

get_r0 <- function(tau, graph, gam=0.5){
  # function for continous time markov sir model with gamma=1
  r0=tau/(tau+gam) * mean(degree(graph)^2 - degree(graph))/mean(degree(graph))
  return(round(r0, 1))
}

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

###### ----------- SIR w/ variable susceptibility 100 REPS-----------  ###### 

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




###### ----------- SIR w/ variable susceptibility (300 REPS) -----------  ###### 

rs <- c(0, 0.3, 0.6, 0.9)
ns <- 1e3
reps <- 1:300
tau <- c(0.02, 0.05, 0.08, 0.16, .24)
alpha <- c(1, 1.5, 2, 2.5)

vars <- expand.grid(r=rs, n=ns, rep=reps, tau=tau, alpha=alpha)

res <- data.frame(size=NA, rep=NA, tau=NA, alph=NA, r0=NA, 
                  degAssort=NA, clustering=NA, pathLen=NA, n_components=NA,
                  r=NA, q=NA, 
                  time_steps=NA, peak_size=NA, 
                  peak_time=NA, tot_inf=NA, 
                  prev_1=NA, prev_2=NA, prev_ratio=NA)

setwd("~/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-rewiring")

for(i in 1:nrow(vars)){

  s=vars[i, "n"]
  r=vars[i, "r"]
  rep=vars[i, "rep"] 
  tau=vars[i, "tau"]
  alph=vars[i, "alpha"]
  
  ## Network data ###
  g <- read.graph(paste0("networks3/G_",
                         r, "N",
                         s, "rep",
                         rep, ".graphml"),
                  format = "graphml")
  res[i, "r0"] <- get_r0(tau, g)
  res[i, "r"] <- r
  res[i, "size"] <- s
  res[i, "rep"] <- rep 
  res[i, "tau"] <- tau
  res[i, "alph"] <- alph
  res[i, "degAssort"] <- assortativity_degree(g, directed = FALSE)
  res[i, "clustering"] <- transitivity(g)
  res[i, "n_components"] <- count_components(g)
  res[i, "pathLen"] <- diameter(g, directed = FALSE)
  res[i, "r_actual"] <- assortativity_nominal(g, types=V(g)$sex, directed=FALSE)
  res[i, "q"] <- modularity(g, membership=V(g)$sex)
  
  if(alph==2) alph="2.0"
  alph=as.character(alph)
  
  ### Epidemic data ### 
  # Combined state information (all time)
  simOut=read.csv(paste0("SIR/var-suscept-3/SIR_R", r, 
                         "_N", s, 
                         "_tau", tau, "_alph", alph,
                         "_rep", rep, ".csv"))
  
  # Final status of nodes 
  nodeOut=read.csv(paste0("SIR/var-suscept-3/Final_R", vars[i, "r"],
                          "_N", s,
                          "_tau", tau, "_alph", alph,
                          "_rep", rep, ".csv"),
                   header = FALSE, col.names=c("node", "state"))
  
  nodeOut$node <- gsub("n", "", nodeOut$node)
  nodeOut$node <- as.numeric(nodeOut$node)
  
  nodeOut <- dplyr::arrange(nodeOut, node)
  #if(nrow(nodeOut) != length(V(g)$sex)) next 
  
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

write.csv(res, "SIR/variable_results_3.csv") # has initial infected =25; reps =300

###### ----------- SIR &  SIRS combined (300 REPS) -----------  ###### 

rs <- c(0, 0.3, 0.6, 0.9)
ns <- 1e3
reps <- 1:300
tau <- c(0.025, 0.045, 0.08, 0.12)
alpha <- c("1.0", "2.0", "3.0")
sig <- c(0, 0.1, 0.2)

vars <- expand.grid(r=rs, n=ns, rep=reps, tau=tau, alpha=alpha, sig=sig)

res <- data.frame(size=rep(NA, nrow(vars)), rep=rep(NA, nrow(vars)), 
                  tau=rep(NA, nrow(vars)), alph=rep(NA, nrow(vars)), 
                  sigma=rep(NA, nrow(vars)), r0=rep(NA, nrow(vars)), 
                  degAssort=rep(NA, nrow(vars)), clustering=rep(NA, nrow(vars)), 
                  pathLen=rep(NA, nrow(vars)), n_components=rep(NA, nrow(vars)),
                  r=rep(NA, nrow(vars)), q=rep(NA, nrow(vars)), 
                  time_steps=rep(NA, nrow(vars)), peak_size=rep(NA, nrow(vars)), 
                  peak_time=rep(NA, nrow(vars)), tot_inf=rep(NA, nrow(vars)), 
                  prev_m=rep(NA, nrow(vars)), prev_f=rep(NA, nrow(vars)), 
                  prev_ratio=rep(NA, nrow(vars)))

setwd("~/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-rewiring")
tmp <- rep(NA, nrow(vars))

for(i in 1:nrow(vars)){
  if (i%%1000==0) print(i)
  s=vars[i, "n"]
  r=vars[i, "r"]
  rep=vars[i, "rep"] 
  tau=vars[i, "tau"]
  alph=vars[i, "alpha"]
  sigm=vars[i, "sig"]
  
  ## Network data ###
  g <- read.graph(paste0("networks3/G_",
                         r, "N",
                         s, "rep",
                         rep, ".graphml"),
                  format = "graphml")
  res[i, "r0"] <- get_r0(tau, g)
  res[i, "r"] <- r
  res[i, "size"] <- s
  res[i, "rep"] <- rep
  res[i, "tau"] <- tau
  res[i, "alph"] <- alph
  res[i, "sigma"] <- sigm
  res[i, "degAssort"] <- assortativity_degree(g, directed = FALSE)
  res[i, "clustering"] <- transitivity(g)
  res[i, "n_components"] <- count_components(g)
  res[i, "pathLen"] <- diameter(g, directed = FALSE)
  res[i, "r_actual"] <- assortativity_nominal(g, types=V(g)$sex, directed=FALSE)
  res[i, "q"] <- modularity(g, membership=V(g)$sex)
  
  alph=as.character(alph)
  
  ### Epidemic data ### 
  # Combined state information (all time)
  simOut=read.csv(paste0("SIRS/SIRS2_R", r, 
                         "_N", s, 
                         "_tau", tau, "_alph", alph, "_sig", sigm, 
                         "_rep", rep, ".csv"))
  
  # Last value for state information (last time point for all ndoes)
  # simulation variables
  res[i, "time_steps"] <- tail(simOut$t, 1) # SIR will die out
  res[i, "peak_size"] <- max(simOut$I.f + simOut$I.m)
  res[i, "peak_time"] <- simOut$t[which.max(simOut$I.f + simOut$I.m)]
  res[i, "tot_inf"] <- ifelse(sigm==0, tail(simOut$R.f + simOut$R.m, 1), mean(simOut$I.f[simOut$t>100] + simOut$I.m[simOut$t>100]))
  tmp[i] <- ifelse(sigm==0, tail(simOut$R.f + simOut$R.m, 1), mean(simOut$I.f[simOut$t>100] + simOut$I.m[simOut$t>100]))
  
  res[i, "prev_m"] <- ifelse(sigm==0, tail(simOut$R.m, 1)/(500), sum(simOut$I.m[simOut$t>100])/(500))
  res[i, "prev_f"] <- ifelse(sigm==0, tail(simOut$R.f, 1)/(500), sum(simOut$I.f[simOut$t>100])/(500))
  res[i, "prev_ratio"] <- res[i, "prev_m"]/res[i, "prev_f"]
}


write.csv(res, "SIRS/results-combined2.csv") # combined SIR and SIRS model results
