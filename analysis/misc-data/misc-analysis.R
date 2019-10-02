library(readxl)
library(tidyverse)
library(ggplot2)
library(GGally)

# Lindsley 2012 COUGH VOLUME BY SEX

dat <- read_xlsx("lindsley-cough-data.xlsx")

dat %>% group_by(sex) %>% summarise(meanVol1=mean(volume1), meanVol2=mean(volume2, na.rm=TRUE))

dat %>% ggplot(aes(volume1)) + geom_histogram() + facet_grid(sex~.)

ggpairs(dat[,c("sex","height", "weight","volume1")],
        mapping=ggplot2::aes(colour = sex))
  

ggpairs(dat[,c("sex","height", "weight","volume2")],
        mapping=ggplot2::aes(colour = sex))

# Lindsley 2010 COUGH VOLUME BY SEX

dat <- read.table("Results_S1.txt", sep="\t", na.strings =c( "n/a", " "))

dat <- dat[,c(3, 11, 25:27)]

dat$meanVol <- rowMeans(dat[,3:5], na.rm = TRUE)

dat %>% group_by(V11) %>% summarise(mean(meanVol))


# Kwok 2018 SEX ASSORTATIVITY IN INDIA
dat <- read.csv("kwok-etal-2018-india.csv")

# doesn't have contact sex
  
