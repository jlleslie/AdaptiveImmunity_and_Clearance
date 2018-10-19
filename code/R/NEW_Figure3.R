#Figure 3
#Analyis of 16S Data from adoptive transfer experiment 
#Question: Does reconsituiton of adaptive immunity in mice alter recovery of the community? 

# Start with clean environment
rm(list=ls())
gc()

# Load dependencies (code coppied from Matt Jenior) 
deps <- c('ggplot2', 'vegan', 'reshape2', 'grid');
for (dep in deps){
  if (dep %in% installed.packages()[,'Package'] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep)

#####Figure 3
#### MDS of D1 communites, after infection before transfer of spleen cells  

#Read in filted and subsampled shared file this file has more data than is used for this analysis 
otu.shared<-read.delim(file="Adaptiveimmuneclear_noD40.42.0.03.subsample.0.03.filter.shared", header = T)
otu.shared$label<-NULL
otu.shared$numOtus<-NULL
#removes label and numOtus columns from the dataframe
#reads in the table of all the samples from the adoptive transfer experiment 
adoptrans.grps<-read.delim(file="2016_RAG_adoptivetransfer.grps.accnos", header = F)
#reads in a file of metadata 
meta.data<-read.delim(file="Adoptivetransfer_metadata.txt", header=T, row.names = 1)

#Pull out the data from the adoptive transfer experiment from the shared file 
otu.shared.adoptrans<-otu.shared[otu.shared$Group %in% adoptrans.grps[,1], ]
row.names(otu.shared.adoptrans)<-otu.shared.adoptrans$Group
otu.shared.adoptrans$Group<-NULL
# question: Was the communites of the mice that cleared differenet before adooptive transfer treatment? 
#add metadata file so you can easily pull out samples based on their metat data
otu_meta<-merge(otu.shared.adoptrans,meta.data, by= 'row.names')
#D1 Dataframe 
D1.otu <- otu_meta[otu_meta$Day == "1", ]
row.names(D1.otu)<-D1.otu$Row.names
D1.otu$Row.names<-NULL
# makes D1 samples only dataframe, remakes the rownames actual rownames rather than a new column 
D1.otu.nometa<-D1.otu[,!colnames(D1.otu)%in% colnames(meta.data)]
#removes metadata columns
D1_nmds <- metaMDS(D1.otu.nometa, distance = "bray",k=2, trymax=100)$points
#runs makes MDS of data based on Bray-Curtis distances (takes into account both eveness and richness )
#clumisly I need to add back metadat to the NMDS file to pull out samples by cage name
D1_nmds_meta<-merge(D1_nmds,meta.data, by= 'row.names')
row.names(D1_nmds_meta)<-D1_nmds_meta$Row.names
D1_nmds_meta$Row.names<-NULL
#reassigns row.names 
#pulling out points for each cage
c143 <-D1_nmds_meta[D1_nmds_meta$Cage =="143", 1:2]
c144 <-D1_nmds_meta[D1_nmds_meta$Cage =="144", 1:2]
c145 <-D1_nmds_meta[D1_nmds_meta$Cage =="145", 1:2]
c146 <-D1_nmds_meta[D1_nmds_meta$Cage =="146", 1:2]
c147 <-D1_nmds_meta[D1_nmds_meta$Cage =="147", 1:2]
c150A <-D1_nmds_meta[D1_nmds_meta$Cage =="150A", 1:2]
c150  <-D1_nmds_meta[D1_nmds_meta$Cage =="150", 1:2]
#make new data frame that I can modify the Cage column in order to calculate the centroid 
D1.cent.nmds<-D1_nmds_meta
D1.cent.nmds$Cage<-c( rep("143",4), rep("other", 17 ))
other<- D1.cent.nmds[D1.cent.nmds$Cage=="other", 1:2]
D1_nmds_meta_centoids <- aggregate(cbind(D1.cent.nmds$MDS1, D1.cent.nmds$MDS2)~ D1.cent.nmds$Cage, data=D1.cent.nmds, mean)

##Making the plot 
#Plot D1 data as an empty plot and then add points for each cage 
#color = treatment group they will end up in #shape= cage 
plot(D1_nmds_meta$MDS1, D1_nmds_meta$MDS2, type = "n",xaxt='n', yaxt='n', cex=0.5, las=1,
     xlab ="MDS axis 1", ylab ="MDS axis 2", xlim = c(-3.1,1.1), ylim=c(-2.7,2.7))
box(which = "plot", lty = "solid", col ="grey80", lwd=5)
axis(side = 2, col="grey80", las=1)
axis(side = 1, col="grey80", las=1)
#abline(v=0,col="grey70",lty=2)
#abline(h=0,col="grey70",lty=2 )
points(c143, pch=21, bg='#f91780', cex=2)
points(c144, pch=22, bg='#fa8c17', cex=2)
points(c145, pch=23, bg='#0095a3', cex=2)
points(c146, pch=24, bg='#f91780', cex=2)
points(c147, pch=13,col='#fa8c17', cex=2)
points(c150, pch=9, col='#f91780', cex=2)
points(c150A, pch=8, col='#0095a3', cex=2)
segments(x0=c143$MDS1, y0=c143$MDS2, x1=D1_nmds_meta_centoids[1,2], y1=D1_nmds_meta_centoids[1,3], col='gray30')
segments(x0=other$MDS1, y0=other$MDS2, x1=D1_nmds_meta_centoids[2,2], y1=D1_nmds_meta_centoids[2,3], col='gray30')
text(0.5,-1.3, labels = c("p = 0.002"), pos=2)
##Statsitics 
D1.otu.nometa$Cage<- sapply(strsplit(row.names(D1.otu.nometa), ".D"), "[", 1)
#Makes a column that is just cage for anosim 

#Global differnces (comparing all cages)
anosim(D1.otu.nometa[,1:476], D1.otu.nometa$Cage, permutations=999, distance='bray')
#D1.otu.nometa[,1:476] is all of the columns with OTUs excluding the cage column 477
#anosim(dat = D1.otu.nometa[, 1:476], grouping = D1.otu.nometa$Cage,      permutations = 999, distance = "bray") 
#Dissimilarity: bray 
#ANOSIM statistic R: 0.4563 
#Significance: 0.001 
#Permutation: free
#Number of permutations: 999

#Testing if cage 143 is different than all others 
D1.otu.nometa$Cage<-c( rep("143",4), rep("other", 17 ))
#rename all other cages other to reduce groups from 7 to 2 
anosim(D1.otu.nometa[,1:476], D1.otu.nometa$Cage, permutations=999, distance='bray')
#Call:
#anosim(dat = D1.otu.nometa[, 1:476], grouping = D1.otu.nometa$Cage,      permutations = 999, distance = "bray") 
#Dissimilarity: bray 
#ANOSIM statistic R: 0.7494 
#Significance:  0.002 
#Permutation: free
#Number of permutations: 999

#Correcting P-values for mutiple comparisons
anosimD1_pvals<-c(0.001,0.002)
round(p.adjust(anosimD1_pvals, method = "BH"),3)
#Corrected P-values
# Anosim all groups: 0.002
#Anosim 143 vs other cages: 0.002
#Report both the corrected P value and the R value as a rule of thumb an R value of 0.7 suggests that the groups are different 

