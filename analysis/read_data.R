# reading in data from epidemics on assorted networks
library(igraph)

# graphs from simulations 

ns <- c(5, 10, 15) * 100
rs <- -9:9 / 10
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
  
  ### Network data ###
  # g <- read.graph(paste0("networks/G_", 
  #                        r, "N", 
  #                        s, "rep", 
  #                        rep, ".graphml"),
  #                 format = "graphml")
  # 
  # res[i, "size"] <- s
  # res[i, "rep"] <- rep
  # res[i, "tau"] <- tau
  # res[i, "degAssort"] <- assortativity_degree(g, directed = FALSE)
  # res[i, "clustering"] <- transitivity(g)
  # res[i, "pathLen"] <- diameter(g, directed = FALSE)
  # res[i, "r"] <- assortativity_nominal(g, types=V(g)$sex, directed=FALSE)
  # res[i, "q"] <- modularity(g, membership=V(g)$sex)

  ### Epidemic data ### 
  # Combined state information (all time)
  simOut=read.csv(paste0("SIR/SIR_R",vars[i, "r"], 
                         "_N", s, 
                         "_tau", tau, 
                         "_rep", vars[i, 'rep'], ".csv"))
  
  simOut$tau <- tau
  simOut$rep <- rep
  simOut$r <- r
  sims[[i]] <- simOut
  
  # Final status of nodes 
  # nodeOut=read.csv(paste0("SIR/Final_R",vars[i, "r"], 
  #                         "_N", s, 
  #                         "_tau", tau, 
  #                         "_rep", vars[i, 'rep'], ".csv"),
  #                  header = FALSE, col.names=c("node", "state"))
  # 
  # nodeOut$node <- gsub("n", "", nodeOut$node)
  # 
  # nodeOut <- dplyr::arrange(nodeOut, node)
  # nodeOut$sex <- V(g)$sex
  # 
  # 
  # # Last value for state information (last time point for all ndoes)
  # # simulation variables
  # res[i, "time_steps"] <- tail(simOut$t, 1)
  # res[i, "peak_size"] <- max(simOut$i)
  # res[i, "peak_time"] <- simOut$t[which.max(simOut$i)]
  # res[i, "tot_inf"] <- table(nodeOut[,2])[["R"]]/s
  # 
  # res[i, "prev_1"] <- length(which(nodeOut$state[nodeOut$sex==1]=="R"))/(s/2)
  # res[i, "prev_2"] <- length(which(nodeOut$state[nodeOut$sex==2]=="R"))/(s/2)
  # res[i, "prev_ratio"] <- res[i, "prev_1"]/res[i, "prev_2"]

  
  
}

#write.csv(res, "results.csv")
write.csv(do.call(bind_rows, sims), "sir-simulations.csv")
