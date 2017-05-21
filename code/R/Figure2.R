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

#Including the following fuction
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

## Pull out CFU at time at D26 (23 days post trainsfer of splenocytes)
colonization.D26<-colonization[colonization$Day == "26", ]
#Remove mouse 143_3 and 144_4, mice that got splenocytes but didn't develope IgG (from both colonization and col D26 data frames)
colonization.no1433.1444<-colonization[colonization$Mouse!="143_3" &colonization$Mouse!="144_4",]
colonization.D26.no1433.1444<-colonization.D26[colonization.D26$Mouse!="143_3" &colonization.D26$Mouse!="144_4",]

#Replace the 0 in the mice that had undectable levles with LOD/squarroot of 2
#the LOD is 100 
fill.in.lod<-100/sqrt(2) #fil.in.lod = 70.71068
colonization.D26.no1433.1444$CFU_g<-replace(colonization.D26.no1433.1444$CFU_g,colonization.D26.no1433.1444$CFU_g==0,fill.in.lod)

#Testing if values from each group are from distint distributions
#Pull out CFU/g  measured in recpient mice that got splenocytes from infected donors
col.spleninfect<-colonization.D26.no1433.1444[colonization.D26.no1433.1444$Treatment_2=="infected_splenocytes",]
col.spleninfect_cfu<-c(col.spleninfect$CFU_g)
#Pull out CFU/g  measured in recpient mice that got splenocytes from uninfected donors 
col.splenmock<-colonization.D26.no1433.1444[colonization.D26.no1433.1444$Treatment_2=="mock_splenocytes",]
col.splenmock_cfu<-c(col.splenmock$CFU_g)
#Pull out CFU/g  measured in recpient mice that got splenocytes from uninfected donors 
col.vehicle<-colonization.D26[colonization.D26$Treatment_2=="vehicle",]
col.vehicle_cfu<-c(col.vehicle$CFU_g)

wilcox.test(col.spleninfect_cfu, col.vehicle_cfu, exact=F)
#data:  col.spleninfect_cfu and col.vehicle_cfu
#W = 29, p-value = 0.5595
wilcox.test(col.spleninfect_cfu, col.splenmock_cfu, exact=F)
#data:  col.spleninfect_cfu and col.splenmock_cfu
#W = 22, p-value = 0.8253
wilcox.test(col.vehicle_cfu, col.splenmock_cfu, exact=F)
#data:  col.vehicle_cfu and col.splenmock_cfu
#W = 20, p-value = 0.4113

#Correcting P-values for mutiple comparisons
col_pvals<-c(0.5595,0.8253,0.4113)
round(p.adjust(col_pvals, method = "BH"),3)
# bengamoni-hersh correction 
# CFU at D26 Infected vs Vehicle : 0.825
# CFU at D26 Infected vs Mock: 0.825
# CFU at D26 Mock vs Vehicle: 0.825



##Plotting D26 Data 
colors<-c("infected_splenocytes"="#f91780", "mock_splenocytes"= "#fa8c17", "vehicle"="#0095a3")
#order the dataset so that it will plot vehicle first on the left
colonization.D26.no1433.1444$Treatment_2<-factor(colonization.D26.no1433.1444$Treatment_2, levels = c("vehicle", "mock_splenocytes", "infected_splenocytes"))

d26.plot<-ggplot(colonization.D26.no1433.1444, aes(x=Treatment_2, y=CFU_g, fill=factor(Treatment_2)))+ 
  geom_jitter(aes(shape=21) , width = 0.2, height = 0.01, size = 4)+
  scale_shape_identity()+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="grey50") +
  scale_fill_manual(values = colors) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)), limits=c(10,4*10^8)) +
  labs(y = expression(paste(" CFU ", "per Gram Feces")))+
  geom_hline(aes(yintercept=100), colour = "gray50", size = 1, linetype=2)
#theme with white background
three.a = d26.plot + 
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
    ,axis.text.x=element_blank()
    ,plot.margin = unit(c(1,1,2,1), "lines")
  )

# Add in the lables for treatment groups outside of the plot
gtext.spleninfect<-textGrob("Splenocytes\n(infected donor)", gp = gpar(fontsize = 10))  
gtext.splenmock<-textGrob("Splenocytes\n(uninfected donor)", gp = gpar(fontsize = 10))  
gtext.vhe<-textGrob("Vehicle",gp = gpar(fontsize = 10)) 
gtext.ns <-textGrob("ns", gp = gpar(fontsize = 13))  


three.a  = three.a  + annotation_custom(gtext.ns, xmin=2, xmax=2, ymin = 8.1, ymax= 8.2) + 
  annotate("segment", x = 1, xend = 3, y = 100000000, yend = 100000000, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=2.5, xmax=2.5, ymin=8.35, ymax=8.6) +  
  annotate("segment", x = 2, xend = 3, y = 200000000, yend = 200000000, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=1.5, xmax=1.5, ymin=8.6, ymax=8.8) +
  annotate("segment", x = 1, xend = 2, y =300000000, yend = 300000000, colour = "black", size = 0.7)+
  annotation_custom(gtext.spleninfect, xmin=3, xmax=3,  ymin=0, ymax=0) +
  annotation_custom(gtext.splenmock, xmin=2, xmax=2,  ymin=0, ymax=0) +
  annotation_custom(gtext.vhe, xmin=1, xmax=1, ymin=0, ymax=0)

g = ggplotGrob(three.a)
g$layout$clip[g$layout$name=="panel"] <- "off"
grid.draw(g)



#Plotting All Cages Over time 

colonization.no1433.1444$CFU_g[colonization.no1433.1444$CFU_g == "0"] <-fill.in.lod
#changes 0 values to LoD/sqt(2) so that it is clear they were not detected 

colonization.no1433.1444.noD1<-colonization.no1433.1444[colonization.no1433.1444$Day != "1", ]
#removes Day 1 data from the data set, the data from this day is comprommised 
#because the lables that were placed on the tubes threw off the weights and they were diluted way too much (and don't have a mechanism to fix it)


cfu.treat<-summaryMED(colonization.no1433.1444.noD1, measurevar="CFU_g", metadata=c("Treatment_2","Day"), na.rm=TRUE)

cfu.cage<-summaryMED(colonization.no1433.1444.noD1, measurevar="CFU_g", metadata=c("Cage","Day"), na.rm=TRUE)
group.cols<-c("infected_splenocytes"="#f91780", "mock_splenocytes"= "#fa8c17", "vehicle"="#0095a3",
              "143"="#f91780", "146"="#f91780" , "150"="#f91780", "144"= "#db9204", "147" ="#db9204", "145"="#006670", "150A"= "#006670")

#Plot data
cfu.plot<-ggplot(cfu.treat, aes(x=Day, y=CFU_g, color=factor(Treatment_2)))+ 
  geom_line(size=2.5)+
  geom_line(data=cfu.cage, aes(x=Day, y=CFU_g, color=factor(Cage)), size= 0.55,linetype = 2) +
  scale_color_manual(values=group.cols)+
  scale_x_continuous(breaks = c(3,6,9,12,15,18,21,24), limits =c(3, 26))
#theme with white background
b = cfu.plot + 
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
    ,axis.text.x=element_text(size=11)
    ,plot.margin = unit(c(1,1,2,1), "lines")
  )
b1 = b + labs(y = " CFU per Gram Feces", x = "Day Post Infection")
b2 = b1+ geom_hline(aes(yintercept=100), colour = "gray50", size = 1, linetype=2)
b3 = b2 + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) 
b3

# Add in the lables for treatment groups outside of the plot
gtext.spleninfect<-textGrob("Splenocytes\n(infected donor)", gp = gpar(fontsize = 8))  
gtext.splenmock<-textGrob("Splenocytes\n(uninfected donor)", gp = gpar(fontsize = 8))  
gtext.vhe<-textGrob("Vehicle",gp = gpar(fontsize = 8)) 

three.b  = b3 + 
   annotate("rect", xmin = 15, xmax = 16, ymin = 100000, ymax = 200000, fill = "#0095a3")+  
   annotation_custom(gtext.vhe, xmin=20, xmax=20,  ymin = 5, ymax =5.3 ) +
   annotate("rect", xmin = 15, xmax = 16, ymin = 10000, ymax = 20000, fill = "#fa8c17")+  
   annotation_custom(gtext.splenmock, xmin = 21, xmax = 21,  ymin = 4, ymax = 4.3) +
   annotate("rect", xmin = 15, xmax = 16, ymin = 1000, ymax = 2000, fill = "#f91780")+  
  annotation_custom(gtext.spleninfect, xmin = 21, xmax = 21,  ymin =3, ymax = 3.3) 

g3 = ggplotGrob(three.b)
g$layout$clip[g$layout$name=="panel"] <- "off"
grid.draw(g3)




#### MDS of D1 communites 
#Read in required  data 
#Read in filted and subsampled shared file this file has more data than is used for this analysis 
otu.shared<-read.delim(file="Adaptiveimmuneclear_noD40.42.0.03.filter.0.03.subsample.shared", header = T)
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


## Figure 4A:  MDS Ordination of D1 Communites 
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
     xlab ="MDS axis 1", ylab ="MDS axis 2", xlim = c(-2.5,1.5), ylim=c(-1.5,1))
box(which = "plot", lty = "solid", col ="grey80", lwd=5)
axis(side = 2, col="grey80", las=1)
axis(side = 1, col="grey80", las=1)
#abline(v=0,col="grey70",lty=2)
#abline(h=0,col="grey70",lty=2 )
points(c143, pch=21, bg='black', cex=2)
points(c144, pch=22, bg='black', cex=2)
points(c145, pch=23, bg='black', cex=2)
points(c146, pch=24, bg='black', cex=2)
points(c147, pch=25, bg='black', cex=2)
points(c150, pch=7, col='black', cex=2)
points(c150A, pch=8, col='black', cex=2)

segments(x0=c143$MDS1, y0=c143$MDS2, x1=D1_nmds_meta_centoids[1,2], y1=D1_nmds_meta_centoids[1,3], col='gray30')
segments(x0=other$MDS1, y0=other$MDS2, x1=D1_nmds_meta_centoids[2,2], y1=D1_nmds_meta_centoids[2,3], col='gray30')
text(1.5, -1, labels = c("Shape: Cage"), pos=2)
text(1.5,-1.3, labels = c("p = 0.004"), pos=2)
text(1.5,-1.5, labels = c("R: 0.7584"), pos=2)
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
#ANOSIM statistic R: 0.7584 
#Significance: 0.004 
#Permutation: free
#Number of permutations: 999

#Correcting P-values for mutiple comparisons
anosimD1_pvals<-c(0.001,0.004)
round(p.adjust(anosimD1_pvals, method = "BH"),3)
#Corrected P-values
# Anosim all groups: 0.002
#Anosim 143 vs other cages: 0.004 
#Report both the corrected P value and the R value as a rule of thumb an R value of 0.7 suggests that the groups are different 

rm(list=ls())
