# Script to generate and save assorted network
library(igraph)
library(tidyverse)
library(magrittr)

setwd("~/Documents/phd/research-projects/miller-tb-assortativity/analysis/network-generation/op-2-simulations")

# function for finding within and between edges
get_e_type <- function(elRow, types){
  ifelse(diff(types[elRow])==0, "w", "b")
}


# net reps
reps <- 100

# max number of trials to make assorted graphs
max.iter <- 2000

# sizes
ns <- c(5e2, 1e3, 1.5e3)
rs <- seq(-9, 9, by=1)/10
alph <- 0.2
eps <- 0.035

vars <- expand.grid(size=ns, assort=rs, alph=alph, eps=eps)

for(v in 1:nrow(vars)){
  print(vars[v, ])
  
  for(k in 1:reps){
    # network params
    s <- vars[v, "size"]
    r_f <- vars[v, "assort"]  # desired assort value
    
    # random network
    Gg <- sample_pa(n=s, # number of vertices
                    power=1, # preferential attach.; linear=1
                    m=5, 
                    directed = FALSE)
    
    males <- sample(1:s, (s/2))
    females <- setdiff(1:s, males)
    sexes <- data.frame(num=c(males, females), sex=c(rep(1, (s/2)), rep(2, (s/2)))) %>% 
      arrange(num)
    
    V(Gg)$sex <- sexes$sex
    
    if(r_f==0) {
      write_graph(Gg,file=paste0("G_", r_f, "N", s, "rep", k, ".graphml"), 
                                                format="graphml")
      next()}
    
    rt <- assortativity_nominal(Gg, types=V(Gg)$sex) # basically 0
    Gg0 <- Gg # copy for plotting
    
    # algorithm params
    alph <- vars[v, "alph"] # re-wire % of between edges at once
    eps <- vars[v, "eps"]  # tolerance level for assortativity
    
    w <- rt     # trials to get to r_f

    for(i in 1:max.iter){
      # Rewire edges for assortativity 
      el <- get.edgelist(Gg)
      eTypes <- apply(el, 1, get_e_type, types=V(Gg)$sex)
      
      if(r_f >=0) er <- which(eTypes=="b")
      if(r_f < 0) er <- which(eTypes=="w")
      
      er <- sample(er, size=round(alph * length(er)))
      GgSub <- subgraph.edges(Gg, eids=er)
      
      GgSub <- get.edgelist(rewire(GgSub, keeping_degseq(niter=50, 
                                                         loops=FALSE)))
      el[er, ] <- GgSub
      
      Gg1 <- graph_from_edgelist(el, directed=FALSE)
      
      # Rewire edges if multiples (1)
      el <- get.edgelist(Gg1) # new edge list
      m <- which(which_multiple(Gg1)) # multiple edge ids
      GgSub <- subgraph.edges(Gg, eids=m)
      GgSub <- get.edgelist(rewire(GgSub, each_edge(prob=1, multiple = FALSE, 
                                                    loops=FALSE)))
      el[m, ] <- GgSub
      Gg1 <- graph_from_edgelist(el, directed=FALSE)
      
      # Calculate assortativity
      V(Gg1)$sex <- sexes$sex
      rt1 <- assortativity_nominal(Gg1, types=V(Gg1)$sex)  
      
      # If new assortativity is closer to rf, keep new graph and its stats
      if(r_f >= 0 & rt1 > (r_f + eps)) next # if too pos. assortative -> skip
      if(r_f < 0 & rt1 < (r_f - eps)) next  # if too neg. assortative -> skip
      
      Gg <- Gg1           # else new graph replaces old
      rt <- rt1           # new assortativity coef replaces old
      w <- c(w, rt)       #  document the change in assortativity
      
      

       if((abs(r_f) - abs(rt)) <= eps) write_graph(simplify(Gg), # remove self edges and mutliple edges
                                                   file=paste0("G_", r_f, "N", s, "rep", k, ".graphml"),
                                                   format="graphml")
    }
  }
}

