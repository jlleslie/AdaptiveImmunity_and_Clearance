###Supplmental Figures for Clearance Paper 

#setworking directory
setwd("~/Desktop/AdaptiveImmunity_and_Clearance/data")

library(ggplot2)
library(grid)
library(scales)
library(vegan)

# Modified from soure of function: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
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

# Supp Figure 1: Colonization in WT donor mice. 
donor.col<-read.delim(file= "Adoptivetransfer_donor_colonization.txt", header = T)
#Replace the 100 in the mice that actually had undectable levles with LOD/squarroot of 2
#the LOD is 100 
fill.in.lod<-100/sqrt(2) #fil.in.lod = 70.71068
donor.col$CFU_g<-replace(donor.col$CFU_g,donor.col$CFU_g== 0,fill.in.lod)

##Determine the Median and IQR for CFU grouped by Cage
donor.col.med<-summaryMED(donor.col, measurevar="CFU_g", metadata=c("Cage","Day"), na.rm=TRUE)
colors.don<-c("1500"="grey", "1502"="black")
donor.plot<-ggplot(donor.col.med, aes(x=Day, y= CFU_g, group=Cage, color=factor(Cage)))+  
  geom_point(size=6)+
  scale_color_manual(values = colors.don) +
  geom_errorbar(aes(ymin=firstquart.25per, ymax=thirdquart.75per), width=1, size=0.9)+
  geom_line(size=0.9) 

#theme with white background
SA1 = donor.plot+ 
  #eliminates background, gridlines and key border
  theme(
    panel.background = element_rect(fill = "white", color = "grey80", size = 2)
    ,panel.grid.major = element_line(color = "gray90", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks= element_line(size = 0.6, colour = "grey90")
    ,axis.ticks.length = unit(0.2, "cm")
    ,legend.background = element_blank ()
    ,axis.text.y=element_text(size=13)
    ,axis.title.y=element_text(size=13)
    ,axis.title.x=element_blank()
    ,axis.text.x=element_text(size=11)
  )
SA1 = SA1 + labs(y = "CFU per Gram Feces")
SA2 = SA1 + geom_hline(aes(yintercept=100), colour = "gray10", size = 0.9, linetype=2)
SA3 = SA2 + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))
SA3

#Supp Figure 2: Total IgG By Cage 
IgG<-read.delim(file="TotalIgG_April_5_2017.txt", header = T)
#replacing cage 150 with 150A because these 
#Plotting the data 
#to clearly show points that were not detected (vs detected at LOD), change values to be below LOD line in the figure 
IgG$Total_IgG<-replace(IgG$Total_IgG, IgG$Total_IgG==0, -500) 
#if you want to change the values again you will have to replace fill.in.lod with -500 etc. 

#assign colors to different treatment groups
colors<-c("infected_splenocytes"="#f91780", "mock_splenocytes"= "#fa8c17", "vehicle"="#0095a3")

##order the dataset so that it will plot vehicle first on the left
IgG$Treatment_2<-factor(IgG$Treatment_2, levels = c("vehicle", "mock_splenocytes", "infected_splenocytes"))
shape_cage<-c("143"= 21, "144"=22, "145" =23, "146"=24, "147" =25, "150" =7 , "150A" =8)
#Make a jitter plot of the data 
igg.plot<-ggplot(IgG, aes(x=Treatment_2, y=Total_IgG, color=factor(Treatment_2), fill=factor(Treatment_2),shape=factor(Cage)))+ 
  scale_shape_manual(values= shape_cage)+
  geom_jitter(width = 0.2, height = 0.01, size= 5)+
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_hline(aes(yintercept=1.56), colour = "grey50", size = 1, linetype=2)

two.B = igg.plot + 
  #eliminates background, gridlines and key border
  theme_bw()
two.B  = two.B + labs(y = "Serum IgG ng/ml")
two.B 



#Supp figure3: MDS of RAG and WT mice co-housed (figure 4B, but analyzed by genotype )
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
text(0.31,-0.16, labels = c("p = 0.085"), pos=2)
text(0.31,-0.13, labels = c("R: 0.1224"), pos=2)


#Stats
anosim_pvalue <- c()
for (i in 1:100){
  anosim_pvalue[i] <- anosim(coho.shared[,1:660], coho.shared$Genotype, permutations=999, distance='bray')$signif
  print(i)
}
print(median(anosim_pvalue))


#Supp figure 4: Colonization of other mice included in Random Forest 

#Colonization WT 2013 Experiment
###Read in the data for all experiments 
cfu<-read.table(file='Colonization_Overtime_630_Allexperiments_copy.txt', header=TRUE)
#Note the data is reported such that if no colonies were seen on the 10^-2 plate, that sample was reported as having 100 CFU/g feces ie the LOD 

##2013 Experiment data, only WT mice
cfu$Cage<-as.factor(cfu$Cage)
cfu.exp13<-cfu[grep("71.",cfu$Cage, value =F),]
#pulls out CFU overtime forthe  2013 experiment only, this is based on the fact that those cages started with 71X numbering  

##Remove D13 Data point (due to issue in plating, samples were plated late in the day after collection)
cfu.exp13.d13<-grep("13",cfu.exp13$Day, value =F)
#pulls out data from D13 
cfu.exp13.NoD13<-cfu.exp13[-c(cfu.exp13.d13),]
#removes D13 data from 2013 data 

##Remove Cage 710 and 711 (they were uninfected)
cfu.710<-grep("710",cfu.exp13.NoD13$Cage, value=F)
cfu.exp13.NoD13<-cfu.exp13.NoD13[-c(cfu.710), ]
#removes cage 710
cfu.711<-grep("711",cfu.exp13.NoD13$Cage, value=F)
#pulls out where cage 711 data points are now you have removed 710 data
cfu.exp13.NoD13<-cfu.exp13.NoD13[-c(cfu.711), ]
cfu.exp13.NoD13$CFU_g <- replace(cfu.exp13.NoD13$CFU_g, cfu.exp13.NoD13$CFU_g==100, 25)
#replaces 100 which is LoD with 25 which is easier to see 
##Determine the Median and IQR for  CFU grouped by Treatment group 
cfu.exp13.NoD13<-summaryMED(cfu.exp13.NoD13, measurevar="CFU_g", metadata=c("Cage","Day"), na.rm=TRUE)

#Plot data
exp2013.plotcfu<-ggplot(cfu.exp13.NoD13, aes(x=Day, y= CFU_g , colour= factor(Cage)))+ 
  geom_errorbar(aes(ymin=firstquart.25per, ymax=thirdquart.75per), width=1, size=0.9)+
  geom_line(size=0.9) +
  geom_point (size=3)
#theme with white background
SB = exp2013.plotcfu+ 
  #eliminates background, gridlines and key border
  theme(
    panel.background = element_rect(fill = "white", color = "grey80", size = 2)
    ,panel.grid.major = element_line(color = "gray90", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks= element_line(size = 0.6, colour = "grey90")
    ,axis.ticks.length = unit(0.2, "cm")
    ,legend.background = element_blank ()
    ,axis.text.y=element_text(size=13)
    ,axis.title.y=element_text(size=13)
    ,axis.title.x=element_blank()
    ,axis.text.x=element_text(size=11)
  )
SB1  = SB + labs(y = expression(paste(Log[10], " CFU ", "per Gram Feces")))
SB2 = SB1+ geom_hline(aes(yintercept=100), colour = "gray10", size = 0.9, linetype=2)
SB3 = SB2 + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))
SB3


#Supp Figure 5: Relative Abundace of OTU 3(Akkermansia)


shared<-read.delim(file="Adaptiveimmuneclear_noD40.42.0.03.subsample.0.03.filter.0.03.pick.shared")
shared$label=NULL
shared$numOtus=NULL
shared$Total.seqs=apply(shared[,2:419], 1, sum)
shared$RelAbund.OTU3= (shared[,4]/shared$Total.seqs*100)

shared.OTU3= shared[ ,c(1,421)]
row.names(shared.OTU3)=shared.OTU3$Group
shared.OTU3$Group = NULL
shared.OTU3$Cage=sapply(strsplit(row.names(shared.OTU3), ".D"), "[", 1)
shared.OTU3$Mouse=sapply(strsplit(row.names(shared.OTU3), "D."), "[", 1)
shared.OTU3$Day=sapply(strsplit(row.names(shared.OTU3), "D"), "[", 2)

#Abudance of OTU3 Akkermansia D21 post infection 
shared.OTU3.D21<-shared.OTU3[shared.OTU3$Day=="21",]
treament.metadata<-read.delim(file="D21.IgGposneg.txt",header = F, row.names = 1)
D21.data<- merge(treament.metadata,shared.OTU3.D21, by='row.names')

plot.D21<-ggplot(D21.data, aes(x=V2, y=RelAbund.OTU3, fill=V2)) +
  geom_boxplot() + theme_bw()
igg.pos<-D21.data[D21.data$V2 =="Splenocytes_IgG_positive",3]
splen.igg.neg<-D21.data[D21.data$V2 =="Splenocytes_IgG_negative",3]
veh<-D21.data[D21.data$V2 =="Vehicle_IgG_negative",3]
plot.D21
wilcox.test(splen.igg.neg,igg.pos)
#data:  splen.igg.neg and igg.pos
#W = 9, p-value = 0.5508
wilcox.test(splen.igg.neg,veh)
#data:  splen.igg.neg and veh
#W = 0, p-value = 0.07143
wilcox.test(igg.pos,veh)
#data:  igg.pos and veh
#W = 12, p-value = 0.0199
#Correcting P-values for mutiple comparisons
recip_pvals<-c(0.008087,0.008087, NA)
round(p.adjust(recip_pvals, method = "BH"),3)



#Abudance of OTU3 Akkermansia before any treatment 
shared.OTU3.preabx<-shared.OTU3[shared.OTU3$Day=="neg12",]
shared.OTU3.preabx$Group<-c(rep("B",1 ), rep("A",1 ), rep("B",4 ), rep("A",1 ), rep("C",2),rep("B",3 ),rep("C",1), rep("B",3 ) )
#since at D-12 none of the mice have been treated with anything, I added random groups A-C to each mouse
#where:
# group A are the mice that will evenutally be Splenocytes_IgG_negative
# group B are the mice that will evenutally be Splenocytes_IgG_positive
# group C are the mice that will evenutally be Vehicle_IgG_negative

plot.Preabx<-ggplot(shared.OTU3.preabx,aes(x=Group, y=RelAbund.OTU3, fill=Group)) +
  geom_boxplot() + scale_y_continuous(limits = c(0,50)) + theme_bw()
plot.Preabx
A<-shared.OTU3.preabx[shared.OTU3.preabx$Group =="A",1]
B<-shared.OTU3.preabx[shared.OTU3.preabx$Group =="B",1]
C<-shared.OTU3.preabx[shared.OTU3.preabx$Group =="C",1]
wilcox.test(A,B)
#data:  A and B
#W = 13, p-value = 0.7669
wilcox.test(A,C)
#data:  A and C
#W = 5, p-value = 0.4
wilcox.test(B,C)
#data:  B and C
#W = 26, p-value = 0.1607

