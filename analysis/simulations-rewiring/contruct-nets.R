# Script to generate and save assorted network
library(igraph)
library(tidyverse)
library(magrittr)

setwd("~/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-rewiring/networks3")

# function for finding within and between edges
get_e_type <- function(elRow, types){
  ifelse(diff(types[elRow])==0, "w", "b")
}

# net reps
reps <- 300

# max number of trials to make assorted graphs
max.iter <- 2000

# sizes
ns <- 1e3
rs <- c(0, 0.3, 0.6, 0.9)
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
    set.seed(floor(runif(1)*1000000))
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
    # calculate current assortativity
    rt <- assortativity_nominal(Gg, types=V(Gg)$sex) # basically 0
    
    # get algorithm params
    alph <- vars[v, "alph"] # re-wire % of between edges at once
    eps <- vars[v, "eps"]  # tolerance level for assortativity
    
    w <- rt     # trials to get to r_f
    
    for(i in 1:max.iter){
      # Rewire edges for assortativity 
      el <- get.edgelist(Gg)
      eTypes <- apply(el, 1, get_e_type, types=V(Gg)$sex)
      
      # only working with positive assortativity now but can do negative with if(r_f < 0) er <- which(eTypes=="w")
      er <- which(eTypes=="b")

      # get edges to rewire from sample from between or within edges (depending on r_f)
      er <- sample(er, size=round(alph * length(er)))
      
      # get the subgraph of edges to rewire
      GgSub <- subgraph.edges(Gg, eids=er)
      
      # rewire edges with keeping degree sequence
      GgSub <- get.edgelist(rewire(GgSub, keeping_degseq(niter=50, 
                                                         loops=FALSE)))
      # insert new edges into graph
      el[er, ] <- GgSub
      
      # make graph object
      Gg1 <- graph_from_edgelist(el, directed=FALSE)
      
      # Make simple graph and replace edges that were not simple
      m <- length(which(which_multiple(Gg1)))
      Gg1 <- igraph::simplify(Gg1)
      Gg1 <- Gg1 + edges(sample(V(Gg1), m*2, replace = FALSE))
      
      # Check for multiple components and if true, reject this rewiring and try another
      if(count_components(Gg1)>1) next()
      
      # Calculate assortativity
      V(Gg1)$sex <- sexes$sex
      rt1 <- assortativity_nominal(Gg1, types=V(Gg1)$sex)  
      
      # If new assortativity is closer to rf, keep new graph and its stats
      if(r_f >= 0 & rt1 > (r_f + eps)) next() # if too pos. assortative -> skip
      if(r_f < 0 & rt1 < (r_f - eps)) next()  # if too neg. assortative -> skip
      
      Gg <- Gg1           # else new graph replaces old
      rt <- rt1           # new assortativity coef replaces old
      w <- c(w, rt)       #  document the change in assortativity
      
      if((abs(r_f) - abs(rt)) <= eps) write_graph(simplify(Gg), # remove self edges and mutliple edges
                                                  file=paste0("G_", r_f, "N", s, "rep", k, ".graphml"),
                                                  format="graphml")
    }
  }
}


