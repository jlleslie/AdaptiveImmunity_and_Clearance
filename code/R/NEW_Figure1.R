###Analysis for Figure 1: 
#Co-housing of RAG and WT mice 
# Experiment Question: Is adaptive immunity required for clearance of C. difficile? 

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

setwd("~/Desktop/AdaptiveImmunity_and_Clearance/data")


#Figure 1A
#Analysis of Co-housing Experiment, community after co-housing but before (abx day -15) RAG vs WT mice 
shared<-read.delim(file="Adaptiveimmuneclear_noD40.42.0.03.subsample.0.03.filter.shared", header=T)
shared$label<-NULL
shared$numOtus<-NULL
meta.data<-read.delim(file="Adaptiveimmuneclear_metadata_noD40.42.tsv", header=T, row.names = 1)
meta.data.coho<-meta.data[meta.data$Year=="2014" & meta.data$Treatment_1=="630",]
meta.data.coho.dneg15<-meta.data.coho[meta.data.coho$Day=="-15" & meta.data.coho$Mouse!="181",]
#Exluded mouse 18-1 because it died during the course of the experiment 
samples<-row.names(meta.data.coho.dneg15)
coho.shared<-shared[shared$Group %in% samples,]
row.names(coho.shared)<-coho.shared$Group
coho.shared$Group=NULL

#Making MDS plot of groups of co-housed mice 
#Generate the MDS 
metaMDS(coho.shared, distance = "bray",k=2, trymax=100)$stress
#Stress =  0.1864309 
Dneg15_nmds <- metaMDS(coho.shared, distance = "bray",k=2, trymax=100)$points
Dneg15_nmds_meta<-merge(Dneg15_nmds,meta.data, by= 'row.names')
row.names(Dneg15_nmds_meta)<-Dneg15_nmds_meta$Row.names
Dneg15_nmds_meta$Row.names<-NULL
#reassigns row.names 
##MDS by genoytpe 
#pulling out points for each co-housing grouping 
pts.RAG<-Dneg15_nmds_meta[Dneg15_nmds_meta$Cage =="16" | Dneg15_nmds_meta$Cage =="18", 1:2]
pts.WT<-Dneg15_nmds_meta[Dneg15_nmds_meta$Cage =="978" |  Dneg15_nmds_meta$Cage =="977", 1:2]
coho.shared$Cage=coho.shared$Genotype
coho.shared$Genotype<-c(rep("RAG1KO",8), rep("WT",9))
Dneg15.cent.nmds<-Dneg15_nmds_meta
RAG<- Dneg15.cent.nmds[Dneg15.cent.nmds$Genotype=="RAG1KO", 1:2]
Dneg15nmdsgeno_centroids <- aggregate(cbind(Dneg15.cent.nmds$MDS1, Dneg15.cent.nmds$MDS2)~ Dneg15.cent.nmds$Genotype, data=Dneg15.cent.nmds, mean)
#Plotting the values 
plot(Dneg15_nmds_meta$MDS1, Dneg15_nmds_meta$MDS2, type = "n",xaxt='n', yaxt='n', cex=0.5, las=1,
     xlab ="MDS axis 1", ylab ="MDS axis 2", xlim = c(-0.27,0.32), ylim=c(-0.25,0.29))
box(which = "plot", lty = "solid", col ="grey80", lwd=5)
axis(side = 2, col="grey80", las=1)
axis(side = 1, col="grey80", las=1)
points(pts.RAG, pch=21, bg='white', cex=2)
points(pts.WT, pch=21, bg='black', cex=2)
segments(x0=pts.WT$MDS1, y0=pts.WT$MDS2, x1=Dneg15nmdsgeno_centroids[2,2], y1=Dneg15nmdsgeno_centroids[2,3], col='gray30')
segments(x0=pts.RAG$MDS1, y0=pts.RAG$MDS2, x1=Dneg15nmdsgeno_centroids[1,2], y1=Dneg15nmdsgeno_centroids[1,3], col='gray30')
points(0.31,0.26, pch=21, bg='white', cex=2)
points(0.31,0.23, pch=21, bg='black', cex=2)
text(0.31, 0.26, labels = c("RAG1KO"), pos=2)
text(0.31, 0.23, labels = c("WT"), pos=2)
text(0.31,-0.16, labels = c("p = 0.087"), pos=2)
#Stats
anosim_pvalue <- c()
for (i in 1:100){
  anosim_pvalue[i] <- anosim(coho.shared[,1:660], coho.shared$Genotype, permutations=999, distance='bray')$signif
  print(i)
}
print(median(anosim_pvalue))
#anosim_pvalue 0.087
#note this value will change a bit due to the nature of ANOSIM 



#Figure 1B
#Colonization overtime in co-housed RAGko and WT mice 

cfutime_data<-read.table(file='Colonization_Overtime_630_Allexperiments_copy.txt', header=TRUE)
#Note the data is reported such that if no colonies were seen on the 10^-2 plate, that sample was reported as having 100 CFU/g feces ie the LOD 

##2014 Exeperiment data, RAG and WT mice
cfutime_data<-cfutime_data[cfutime_data$Experiment == "2014",]

#Pulls out Mice that were infected 
cfutime_data.infected<-cfutime_data[cfutime_data$Treatment_1 =="630",]
#Replace the 100 in the mice that actually had undectable levles with LOD/squarroot of 2
#the LOD is 100 
fill.in.lod<-100/sqrt(2) #fil.in.lod = 70.71068
cfutime_data.infected$CFU_g<-replace(cfutime_data.infected$CFU_g,cfutime_data.infected$CFU_g==100,fill.in.lod)


## Modified from soure of function: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
## Gives count, first quartile, median and thrid quartile 
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   metadatas: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's

summaryMED<-function(data=NULL, measurevar, metadata=NULL, na.rm=FALSE, .drop=TRUE){
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  data1<- ddply(data, metadata, .drop=.drop,
                .fun=function(xx, col){
                  c(N  = length2(xx[[col]], na.rm=na.rm),
                    firstquart = (quantile(xx[[col]], na.rm=na.rm)[2]),
                    median =  median (xx[[col]], na.rm=na.rm),
                    thridquart = (quantile(xx[[col]], na.rm=na.rm)[4])
                  )
                }, 
                measurevar
  )
  data1 <- rename(data1, c("median" = measurevar))
  data1 <- rename(data1, c("firstquart.25%" = "firstquart.25per"))
  data1 <- rename(data1, c("thridquart.75%" = "thirdquart.75per"))
  return(data1)
}


##Determine the Median and IQR for CFU grouped by Cage
cfutime_data.infected.med<-summaryMED(cfutime_data.infected, measurevar="CFU_g", metadata=c("Cage","Day"), na.rm=TRUE)
cfutime_data.infected.med$Genotype<-c(rep("RAG",16), rep("WT",16))
#add genotype column 


col=c("16"="#61bfee", "978"="#61bfee","18"="#bb5fa1","977" ="#bb5fa1")
shape_A=as.numeric(c("RAG"="1", "WT" = "19"))
plot4.b<-ggplot(cfutime_data.infected.med, aes(x=Day, y= CFU_g, group=Cage, color=factor(Cage), shape=Genotype))+  
      geom_point(size=6)+
      scale_shape_manual(values=shape_A) +
      geom_errorbar(aes(ymin=firstquart.25per, ymax=thirdquart.75per), width=1, size=0.9)+
      geom_line(size=0.9) +
      scale_color_manual(values=col)
#theme with white background
b = plot4.b+ 
  #eliminates background, gridlines and key border
  theme(
    panel.background = element_rect(fill = "white", color = "grey80", size = 2)
    ,panel.grid.major = element_line(color = "gray90", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks= element_line(size = 0.6, colour = "grey90")
    ,axis.ticks.length = unit(0.2, "cm")
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="none"  #if using this as a single figure change "none" to "top" or "bottom" and remove comment from the following 2 lines
    ,axis.text.y=element_text(size=13)
    ,axis.title.y=element_text(size=13)
    ,axis.title.x=element_blank()
    ,axis.text.x=element_text(size=11)
  )
b1 = b+ labs(y = " CFU per Gram Feces")
b2 = b1+ geom_hline(aes(yintercept=100), colour = "gray10", size = 0.9, linetype=2)
b3 = b2 + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))
b3

cfutime_data.infected$CFU_g<-replace(cfutime_data.infected$CFU_g,cfutime_data.infected$CFU_g==100,fill.in.lod)
#replace placed holder LOD values with LOD/squareroot(2) 
d40.f<-cfutime_data.infected[cfutime_data.infected$Day == "40" & cfutime_data.infected$Gender=="female",11]
d40.m<-cfutime_data.infected[cfutime_data.infected$Day == "40" & cfutime_data.infected$Gender=="male",11]
#note the group that cleared were male and while the female groups remained colonized
wilcox.test(d40.f, d40.m)
#data:  d40.f and d40.m
#W = 72, p-value = 0.0002166



#Figure 1C
#Analysis of Co-housing Experiment, before abx by co-housing group not gneotype 
shared<-read.delim(file="Adaptiveimmuneclear_noD40.42.0.03.subsample.0.03.filter.shared", header=T)
shared$label<-NULL
shared$numOtus<-NULL
meta.data<-read.delim(file="Adaptiveimmuneclear_metadata_noD40.42.tsv", header=T, row.names = 1)
meta.data.coho<-meta.data[meta.data$Year=="2014" & meta.data$Treatment_1=="630",]
meta.data.coho.dneg15<-meta.data.coho[meta.data.coho$Day=="-15" & meta.data.coho$Mouse!="181",]
#Exluded mouse 18-1 because it died during the course of the experiment 
samples<-row.names(meta.data.coho.dneg15)
coho.shared<-shared[shared$Group %in% samples,]
row.names(coho.shared)<-coho.shared$Group
coho.shared$Group=NULL

#Making MDS plot of groups of co-housed mice 
#Generate the MDS 
metaMDS(coho.shared, distance = "bray",k=2, trymax=100)$stress
#Stress =  0.1864309 
Dneg15_nmds <- metaMDS(coho.shared, distance = "bray",k=2, trymax=100)$points
Dneg15_nmds_meta<-merge(Dneg15_nmds,meta.data, by= 'row.names')
row.names(Dneg15_nmds_meta)<-Dneg15_nmds_meta$Row.names
Dneg15_nmds_meta$Row.names<-NULL
#reassigns row.names 

#pulling out points for each co-housing grouping 
pts.A<-Dneg15_nmds_meta[Dneg15_nmds_meta$Cage =="16" | Dneg15_nmds_meta$Cage =="978", 1:2]
pts.B<-Dneg15_nmds_meta[Dneg15_nmds_meta$Cage =="18" |  Dneg15_nmds_meta$Cage =="977", 1:2]
Dneg15.cent.nmds<-Dneg15_nmds_meta
Dneg15.cent.nmds$Cage<-c(rep("A",4), rep("B",8), rep("A",5))
A<- Dneg15.cent.nmds[Dneg15.cent.nmds$Cage=="A", 1:2]
Dneg15nmds_centroids <- aggregate(cbind(Dneg15.cent.nmds$MDS1, Dneg15.cent.nmds$MDS2)~ Dneg15.cent.nmds$Cage, data=Dneg15.cent.nmds, mean)

#Plotting the values 
plot(Dneg15_nmds_meta$MDS1, Dneg15_nmds_meta$MDS2, type = "n",xaxt='n', yaxt='n', cex=0.5, las=1,
     xlab ="MDS axis 1", ylab ="MDS axis 2", xlim = c(-0.27,0.32), ylim=c(-0.25,0.29))
box(which = "plot", lty = "solid", col ="grey80", lwd=5)
axis(side = 2, col="grey80", las=1)
axis(side = 1, col="grey80", las=1)
points(pts.A, pch=21, bg='#61bfee', cex=2)
points(pts.B, pch=21, bg='#bb5fa1', cex=2)
segments(x0=pts.B$MDS1, y0=pts.B$MDS2, x1=Dneg15nmds_centroids[2,2], y1=Dneg15nmds_centroids[2,3], col='gray30')
segments(x0=pts.A$MDS1, y0=pts.A$MDS2, x1=Dneg15nmds_centroids[1,2], y1=Dneg15nmds_centroids[1,3], col='gray30')
points(0.31,0.26, pch=21, bg='#61bfee', cex=2)
points(0.31,0.23, pch=21, bg='#bb5fa1', cex=2)
text(0.31, 0.26, labels = c("Group A"), pos=2)
text(0.31, 0.23, labels = c("Group B"), pos=2)
text(0.31,-0.26, labels = c("p = 0.047"), pos=2)
#for the ANOSIM you need to assign samples to either group A (cage 16 and 978) or B (18 and 977) 
coho.shared$Cage<-sapply(strsplit(row.names(coho.shared), ".D"), "[", 1)
coho.shared$Cage<-c(rep("A",4), rep("B",8), rep("A",5))
#Assigns cage 16 and 978 to group A and 16 and 977 to group B 

#Loop runs anosim itteratively 100 times and reports the median value 
anosim_pvalue <- c()
for (i in 1:100){
  anosim_pvalue[i] <- anosim(coho.shared[,1:660], coho.shared$Cage, permutations=999, distance='bray')$signif
  print(i)
}
print(median(anosim_pvalue))
#Significance:  0.047
#Note that this value will change a bit every time it is run 






