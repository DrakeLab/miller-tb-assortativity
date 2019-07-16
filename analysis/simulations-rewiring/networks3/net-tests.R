# analyze networks created by rewiring algorithms
library(igraph)
library(tidyverse)

##### --- 1. networks without correctings for components (possibly not mulitple edges) ####
ns <- c(5, 10, 15) * 100
rs <- 0:9 / 10
reps <- 1:100
vars <- expand.grid(r=rs, n=ns, rep=reps)

res1 <- data.frame(rep=NA, degAssort=NA, clustering=NA, mdeg=NA,
                  pathLen=NA, nComponents=NA, simple=NA, giant=NA, nodes=NA,edges=NA, 
                  r=NA, q=NA)

setwd("~/Documents/phd/res1earch-projects/miller-tb-assortativity/analysis/simulations-rewiring")

for(i in 1:nrow(vars)){
  s=vars[i, "n"]
  r=vars[i, "r"]
  rep=vars[i, "rep"]

  ## Network data ###
  g <- read.graph(paste0("networks/G_",
                         r, "N",
                         s, "rep",
                         rep, ".graphml"),
                  format = "graphml")
  
  res1[i, "rep"] <- rep
  res1[i, "nodes"] <- vcount(g)
  res1[i, "edges"] <- ecount(g)
  res1[i, "mdeg"] <- mean(degree(g))
  res1[i, "simple"] <- is_simple(g)
  res1[i, "nComponents"] <- count_components(g)
  res1[i, "giant"] <- max(components(g)$csize, na.rm=TRUE)
  res1[i, "degAssort"] <- assortativity_degree(g, directed = FALSE)
  res1[i, "clustering"] <- transitivity(g)
  res1[i, "pathLen"] <- diameter(g, directed = FALSE)
  res1[i, "r"] <- assortativity_nominal(g, types=V(g)$sex, directed=FALSE)
  res1[i, "q"] <- modularity(g, membership=V(g)$sex)
}
res1$version=1
#saveRDS(res1, file='v1nets.rds')


##### --- 2. networks with correcting for components #####

# networks with correcting for components
rs <- c(0, 0.3, 0.6, 0.9)
ns <- 1e3
reps <- 1:300

vars <- expand.grid(r=rs, n=ns, rep=reps)

res2 <- data.frame(rep=NA, degAssort=NA, clustering=NA, mdeg=NA,
                          pathLen=NA, nComponents=NA, simple=NA, giant=NA, nodes=NA,edges=NA, 
                          r=NA, q=NA)

setwd("~/Documents/phd/research-projects/miller-tb-assortativity/analysis/simulations-rewiring")

for(i in 1:nrow(vars)){
  s=vars[i, "n"]
  r=vars[i, "r"]
  rep=vars[i, "rep"]
  tau=vars[i, "tau"]
  alph=vars[i, "alpha"]
  
  ## Network data ###
  g <- read.graph(paste0("networks2/G_",
                         r, "N",
                         s, "rep",
                         rep, ".graphml"),
                  format = "graphml")
  
  res2[i, "rep"] <- rep
  res2[i, "nodes"] <- vcount(g)
  res2[i, "edges"] <- ecount(g)
  res2[i, "mdeg"] <- mean(degree(g))
  res2[i, "simple"] <- is_simple(g)
  res2[i, "nComponents"] <- count_components(g)
  res2[i, "giant"] <- max(components(g)$csize, na.rm=TRUE)
  res2[i, "degAssort"] <- assortativity_degree(g, directed = FALSE)
  res2[i, "clustering"] <- transitivity(g)
  res2[i, "pathLen"] <- diameter(g, directed = FALSE)
  res2[i, "r"] <- assortativity_nominal(g, types=V(g)$sex, directed=FALSE)
  res2[i, "q"] <- modularity(g, membership=V(g)$sex)
}
res2$version=2
#saveRDS(res2, file='v2nets.rds')


##### --- Comparison of algorithms 
compareNets <- bind_rows(readRDS("v2nets.rds"), readRDS("v1nets.rds"))%>%
  mutate(roundR=round(r, 1))

compareNets %>%
  filter(nodes==1000) %>%
  select(-simple, -nodes, -edges) %>%
  gather(stat, val, 2:6) %>%
  ggplot(aes(factor(roundR), val)) + 
  geom_boxplot() + facet_grid(stat~ version, scales="free_y")

compareNets %>%
  filter(version==2) %>%
  gather(stat, val, c(2:10, 13)) %>%
  ggplot(aes(factor(roundR), val)) + 
  geom_boxplot() + facet_wrap(stat~ ., scales="free_y")
