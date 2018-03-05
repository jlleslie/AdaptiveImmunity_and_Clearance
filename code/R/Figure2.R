#Figure 2: Colonization over time  in adapotive transfer mice and D1 MDS 
#For this figure and all main text figures, the two mice that did not develop IgG were excluded from the analysis. 

# Start with clean environment
rm(list=ls())
gc()

setwd("~/Desktop/AdaptiveImmunity_and_Clearance/data")
#read in the data 
colonization<-read.table(file="Colonization_RAG_Recips.txt", header=T)


# Load dependencies (code coppied from Matt Jenior) 
deps <- c('ggplot2', 'vegan', 'reshape2', 'grid', 'scales');
for (dep in deps){
  if (dep %in% installed.packages()[,'Package'] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep)

#


