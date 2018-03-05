#Figure 3
#Analyis of 16S Data from adaoptive Transfer Experiment 
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

#####Figure 3A 
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

## Figure 3B: Bray-Curtis Dissimparity between Before Abx to D21 post Infection
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

#add metadata file so you can easily pull out samples based on their metat data
#remove mice 1434 and 1444 because they didn't have sucessful transfer of adaptive immunity 
otu_meta<-merge(otu.shared.adoptrans,meta.data, by= 'row.names')
otu_meta.no1443.1444<-otu_meta[otu_meta$Mouse!="1433" & otu_meta$Mouse!="1444",]
# Preabx and D21 Dataframe 
Dneg12.D21.otu<-otu_meta.no1443.1444[otu_meta.no1443.1444$Day == "-12" | otu_meta.no1443.1444$Day =="21", ]
#note otu_meta was filtered 
row.names(Dneg12.D21.otu)<-Dneg12.D21.otu$Row.names
Dneg12.D21.otu$Row.names<-NULL
#makes Dneg12 and D21 only data frame
#when you merge the row.names becomes a column itself, the above 2 lines returns it to be actual rownames
Dneg12.D21.otu.nometa<-Dneg12.D21.otu[,!colnames(Dneg12.D21.otu)%in% colnames(meta.data)]
#remove the metadata

#calulating distance matrix between all Dneg12 and D21 samples based on bray-curtis metric
Dneg12.D21.dist<-vegdist(Dneg12.D21.otu.nometa, method = "bray", upper=T, diag=T)
Dneg12.D21.dist<-as.matrix(Dneg12.D21.dist)

#Make a vector of all mice in each treatment group exluding the mice  didn't have detectable IgG at D26pi (see figure 2) 
spl_infect<-c("1432" ,"1434", "1461","1503" ,"1504" ,"1505")
    #Don't include  mouse 1433
spl_mock<-c("1441", "1442", "1443" , "1471", "1472")
  #Don't includee mouse 1444
veh<-c( "1452" ,"1454" ,"1501")

# This fuction requres a vector of characters (the mice) and the data.frame. It finds all the values that compare D1 to D21 for each mouse
get.data.Dneg12<-function(data, vector_mice){
  output = c()
  output_index = 1
  for (row_str in row.names(data)) {
    row_vect = unlist(strsplit(row_str,"D"))
    if (row_vect[2] != "neg12"){
      next 
    }
    else if (row_vect[1] %in% vector_mice){
      col_str = paste(row_vect[1],"D21", sep = "")
      output[output_index] = data[row_str, col_str]
      output_index = output_index + 1
    }
  } 
  return(output)  
}

#Running the get.data.Dneg12 function 
dist.infect = get.data.Dneg12(Dneg12.D21.dist, spl_infect)
dist.mock = get.data.Dneg12(Dneg12.D21.dist, spl_mock)
dist.veh =get.data.Dneg12(Dneg12.D21.dist, veh)

#Testing if the Bray-Curtis Dissimiarities are different between the treatment groups. 
wilcox.test(dist.infect, dist.mock)
#data:  dist.infect and dist.mock
#W = 5, p-value = 0.08225
wilcox.test(dist.veh, dist.mock)
#data:  dist.veh and dist.mock
#W = 2, p-value = 0.1429
wilcox.test(dist.veh, dist.infect)
#data:  dist.veh and dist.infect
#W = 7, p-value = 0.7143

#Correcting P-values for mutiple comparisons
recip_pvals<-c(0.08225, 0.1429,0.7143)
round(p.adjust(recip_pvals, method = "BH"),3)
# D1 to D21 disimilarity in mice getting spleen cells from mock vs infected animals:  0.214
# D1 to D21 disimilarity in mice getting spleen cells from mock infected animals vs vehicle: 0.214
# D1 to D21 disimilarity in mice getting spleen cells from infected animals vs vehicle: 0.714

##### Testing if there is a difference if you remove all of the mice that were in cage 143(started out at D1 different than others)
spl_infect.n143<-c("1461","1503" ,"1504" ,"1505")
#Don't include  mouse 1433
dist.infect.n143 = get.data.Dneg12(Dneg12.D21.dist, spl_infect.n143)
wilcox.test(dist.infect.n143, dist.mock)
#W = 5, p-value = 0.2857
wilcox.test(dist.veh, dist.infect.n143)
#W = 1, p-value = 0.1143
#conclusion: removing those cages doesn't alter the fact that there is no significant difference between the groups 

#Plotting the results including cage 143 (cage that cleared)

#Adding NA to make a table to graph 
dist.veh.1 = c(dist.veh, NA, NA,NA)
dist.mock.1 = c(dist.mock, NA)
dneg12vd21.dist<-as.data.frame(cbind(dist.veh.1, dist.mock.1, dist.infect))

#rquires package reshape2 
Dneg12.D21dists.long<-melt(dneg12vd21.dist)
#assign colors to different treatment groups
colors<-c("dist.infect"="#f91780", "dist.mock.1"= "#fa8c17", "dist.veh.1"="#0095a3")

bc.plot<-ggplot(Dneg12.D21dists.long, aes(x=variable, y=value, fill=factor(variable)))+ 
  geom_jitter(aes(shape=21) , width = 0.2, height = 0.01, size = 5)+
  scale_shape_identity()+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="grey50") +
  scale_fill_manual(values = colors) +
  scale_y_continuous( limits = c(0, 1.1)) 

fig3B = bc.plot + 
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
fig3B = fig3B + labs(y = "Bray – Curtis dissimilarity\nPre-Antibiotics to Day 21 post infection")
fig3B

# Add in the lables for treatment groups outside of the plot
gtext.spleninfect<-textGrob("Splenocytes\n(infected donor)", gp = gpar(fontsize = 10))  
gtext.splenmock<-textGrob("Splenocytes\n(uninfected donor)", gp = gpar(fontsize = 10))  
gtext.vhe<-textGrob("Vehicle",gp = gpar(fontsize = 10)) 
gtext.ns<-textGrob("ns", gp = gpar(fontsize = 13)) 
fig3B   = fig3B   + annotation_custom(gtext.ns, xmin=1.5, xmax=1.5, ymin=1.02, ymax =1.02) +  #adding ns for comparsion between splenocytes-infect vs vehicle 
  annotate("segment", x = 1, xend = 2, y = 1, yend = 1, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=2, xmax=2, ymin=1.06, ymax=1.06) +  #adding ns for comparsion between splenocytes- uninfect vs vehicle 
  annotate("segment", x = 1, xend = 3, y =1.04, yend = 1.04, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=2.5, xmax=2.5, ymin=1.1, ymax=1.1) + #adding ns for comparsion between splenocytes- uninfect vs plenocytes- uninfect
  annotate("segment", x = 2, xend =3, y = 1.08, yend = 1.08, colour = "black", size = 0.7)+
  annotation_custom(gtext.spleninfect, xmin=3, xmax=3,  ymin=-0.09, ymax=-0.08) +
  annotation_custom(gtext.splenmock, xmin=2, xmax=2, ymin=-0.09, ymax=-0.08) +
  annotation_custom(gtext.vhe, xmin=1, xmax=1, ymin=-0.09, ymax=-0.08)
g = ggplotGrob(fig3B)
g$layout$clip[g$layout$name=="panel"] <- "off"
grid.draw(g)



# Figure 3C: D21 Comparing inverse simpson divesity between treatment croups 

D21.shared<-read.delim(file="D21_adoptivetrans.no1433.1444.0.03.subsample.0.03.filter.0.03.pick.shared")
#reads in the full D21 sample shared file
row.names(D21.shared) = D21.shared$Group
D21.shared$label <- NULL
D21.shared$Group <- NULL
D21.shared$numOtus <- NULL
D21.shared.div<-as.matrix(diversity(D21.shared, index = "invsimpson"))
meta.data.D21<-meta.data[meta.data$Day =="21",]
D21.invsimp<-as.data.frame(merge(meta.data.D21, D21.shared.div, by = "row.names"))
D21.invsimp.mod<-D21.invsimp[D21.invsimp$Mouse!="1433" & D21.invsimp$Mouse!="1444", ]
#removes the mice who didn't develope IgG 
colors<-c("infected_splenocytes"="#f91780", "mock_splenocytes"= "#fa8c17", "vehicle"="#0095a3")
D21.invsimp.mod$Treatment_2<-factor(D21.invsimp.mod$Treatment_2, levels = c("vehicle", "mock_splenocytes", "infected_splenocytes"))

#Stats
infect.d21<-D21.invsimp.mod[D21.invsimp.mod$Treatment_2=="infected_splenocytes",12]
mock.d21<-D21.invsimp.mod[D21.invsimp.mod$Treatment_2=="mock_splenocytes",12]
veh.d21<-D21.invsimp.mod[D21.invsimp.mod$Treatment_2=="vehicle",12]

wilcox.test(infect.d21, mock.d21)
#data:  infect.d21 and mock.d21
#W = 21, p-value = 0.9433
wilcox.test(infect.d21, veh.d21)
#data:  infect.d21 and veh.d21
#W = 27, p-value = 0.7546
wilcox.test(mock.d21, veh.d21)
#data:  mock.d21 and veh.d21
#W = 18, p-value = 0.6623
#Correction for mutiple comparisons 
three.b.pval<-c(0.9433, 0.7546, 0.6623)
round(p.adjust(three.b.pval, method = "BH"),3)
#[1]  0.943 0.943 0.943
#Plotting Inv Simpson Divserity D21 removing mice that got sleenocytes but didn't have IgG
plot<-ggplot(D21.invsimp.mod, aes(x=Treatment_2, y=V1, fill=factor(Treatment_2) ))+ 
  geom_jitter(aes(shape=21) , width = 0.2, height = 0.01, size = 4)+
  scale_shape_identity()+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="grey50") +
  scale_fill_manual(values = colors) +
  labs(y = "Inverse Simpson Index")+
   scale_y_continuous( limits = c(0,25)) 
  
fig3C= plot + 
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
fig3C
gtext.spleninfect<-textGrob("Splenocytes\n(infected donor)", gp = gpar(fontsize = 10))  
gtext.splenmock<-textGrob("Splenocytes\n(uninfected donor)", gp = gpar(fontsize = 10))  
gtext.vhe<-textGrob("Vehicle",gp = gpar(fontsize = 10)) 
gtext.ns <-textGrob("ns", gp = gpar(fontsize = 13)) 
fig3C  = fig3C   + annotation_custom(gtext.ns, xmin=1.5, xmax=1.5, ymin=21, ymax =21) +  #adding ns for comparsion between splenocytes-infect vs vehicle 
  annotate("segment", x = 1, xend = 2, y = 20, yend = 20, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=2, xmax=2, ymin=23, ymax=23) +  #adding ns for comparsion between splenocytes- uninfect vs vehicle 
  annotate("segment", x = 1, xend = 3, y =22, yend = 22, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=2.5, xmax=2.5, ymin=25, ymax=25) + #adding ns for comparsion between splenocytes- uninfect vs plenocytes- uninfect
  annotate("segment", x = 2, xend =3, y = 24, yend = 24, colour = "black", size = 0.7)+
  annotation_custom(gtext.spleninfect, xmin=3, xmax=3,  ymin=-3.5, ymax=-3.5) +
  annotation_custom(gtext.splenmock, xmin=2, xmax=2, ymin=-3.5, ymax=-3.5) +
  annotation_custom(gtext.vhe, xmin=1, xmax=1, ymin=-3.5, ymax=-3.5)
g.3c = ggplotGrob(fig3C)
g.3c$layout$clip[g.3c$layout$name=="panel"] <- "off"
grid.draw(g.3c)



