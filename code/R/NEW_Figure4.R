###Analysis for Figure 4: 
#Analysis of the gut bacterial community of reconstituted RAG1 mice 

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

setwd("~/Desktop/repos/AdaptiveImmunity_and_Clearance/data")
## Figure 4A: Bray-Curtis Dissimparity between Before Abx to D21 post Infection
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

fig4A = bc.plot + 
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
fig4A = fig4A+ labs(y = "Bray – Curtis dissimilarity\nPre-Antibiotics to Day 21 post infection")
fig4A

# Add in the lables for treatment groups outside of the plot
gtext.spleninfect<-textGrob("Splenocytes\n(infected donor)", gp = gpar(fontsize = 10))  
gtext.splenmock<-textGrob("Splenocytes\n(uninfected donor)", gp = gpar(fontsize = 10))  
gtext.vhe<-textGrob("Vehicle",gp = gpar(fontsize = 10)) 
gtext.ns<-textGrob("ns", gp = gpar(fontsize = 13)) 
fig4A  = fig4A   + annotation_custom(gtext.ns, xmin=1.5, xmax=1.5, ymin=1.02, ymax =1.02) +  #adding ns for comparsion between splenocytes-infect vs vehicle 
  annotate("segment", x = 1, xend = 2, y = 1, yend = 1, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=2, xmax=2, ymin=1.06, ymax=1.06) +  #adding ns for comparsion between splenocytes- uninfect vs vehicle 
  annotate("segment", x = 1, xend = 3, y =1.04, yend = 1.04, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=2.5, xmax=2.5, ymin=1.1, ymax=1.1) + #adding ns for comparsion between splenocytes- uninfect vs plenocytes- uninfect
  annotate("segment", x = 2, xend =3, y = 1.08, yend = 1.08, colour = "black", size = 0.7)+
  annotation_custom(gtext.spleninfect, xmin=3, xmax=3,  ymin=-0.09, ymax=-0.08) +
  annotation_custom(gtext.splenmock, xmin=2, xmax=2, ymin=-0.09, ymax=-0.08) +
  annotation_custom(gtext.vhe, xmin=1, xmax=1, ymin=-0.09, ymax=-0.08)
g = ggplotGrob(fig4A)
g$layout$clip[g$layout$name=="panel"] <- "off"
grid.draw(g)



# Figure 4B: D21 Comparing inverse simpson divesity between treatment croups 

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

fig4B= plot + 
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
fig4B
gtext.spleninfect<-textGrob("Splenocytes\n(infected donor)", gp = gpar(fontsize = 10))  
gtext.splenmock<-textGrob("Splenocytes\n(uninfected donor)", gp = gpar(fontsize = 10))  
gtext.vhe<-textGrob("Vehicle",gp = gpar(fontsize = 10)) 
gtext.ns <-textGrob("ns", gp = gpar(fontsize = 13)) 
fig4B  = fig4B   + annotation_custom(gtext.ns, xmin=1.5, xmax=1.5, ymin=21, ymax =21) +  #adding ns for comparsion between splenocytes-infect vs vehicle 
  annotate("segment", x = 1, xend = 2, y = 20, yend = 20, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=2, xmax=2, ymin=23, ymax=23) +  #adding ns for comparsion between splenocytes- uninfect vs vehicle 
  annotate("segment", x = 1, xend = 3, y =22, yend = 22, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=2.5, xmax=2.5, ymin=25, ymax=25) + #adding ns for comparsion between splenocytes- uninfect vs plenocytes- uninfect
  annotate("segment", x = 2, xend =3, y = 24, yend = 24, colour = "black", size = 0.7)+
  annotation_custom(gtext.spleninfect, xmin=3, xmax=3,  ymin=-3.5, ymax=-3.5) +
  annotation_custom(gtext.splenmock, xmin=2, xmax=2, ymin=-3.5, ymax=-3.5) +
  annotation_custom(gtext.vhe, xmin=1, xmax=1, ymin=-3.5, ymax=-3.5)
g.3c = ggplotGrob(fig4B)
g.3c$layout$clip[g.3c$layout$name=="panel"] <- "off"
grid.draw(g.4B)





##Figure 4C Finding OTUS that discriminate IgG+ vs IgG- (vehicle) mice
#For this analysis, using mothur version 1.39.0, from the shared file "Adaptiveimmuneclear_noD40.42.0.03.subsample.0.03.filter.shared,"
#The shared file had already been subsampled to 10,000 sequences and then filtered so that each OTU was in at least 6 samples (lowest group n)
#I pulled out D21 samples the shared file except for samples 1433D21 and 1444D21 because they didn't have IgG despite getting splenocytes
# I renamed this new shared file "D21_adoptivetrans.no1433.1444.0.03.subsample.0.03.filter.0.03.pick.shared"
# I ran lefse using this new shared file using a  design file had sample  and if they were IgG + or not

######Plot Using Log Transformed Data NOT Relative Abundace ####################
#Using lefse file without Samples from 1443 and 1444 because they are unclear phenotype interms of adaptive Immunity
#Read in lefse result
lefse<-read.delim(file="D21_adoptivetrans.no1433.1444.0.03.subsample.0.03.filter.0.03.pick.0.03.lefse_summary", header=T)
lefse.otus.df<- na.omit(lefse[lefse$LDA >="2",])
#removes all the OTUS associaed with empty values and LDA<2
lefse.otus.df<-lefse.otus.df[order(-lefse.otus.df$LDA), ]
#orders the dataframe by LDA value
lefse.otus<-lefse.otus.df[1:10,]
D21.shared<-read.delim(file="D21_adoptivetrans.no1433.1444.0.03.subsample.0.03.filter.0.03.pick.shared")
#reads in the full D21 sample shared file
row.names(D21.shared) = D21.shared$Group
D21.shared$label <- NULL
D21.shared$Group <- NULL
D21.shared$numOtus <- NULL

# Function to subsample shared files
rarefyOTU <- function(shared, subSize) {
  shared <- t(shared)
  for (x in 1:ncol(shared)) {
    shared[,x] <- as.vector(rrarefy(shared[,x], sample=subSize))
  }
  shared <- as.data.frame(t(shared))
  return(shared)
}

#D21.shared.sub <- rarefyOTU(D21.shared, 5000)
#D21.shared.log<- log10(D21.shared.sub + 1)
D21.shared.log<- log10(D21.shared + 1) #calcs Log10 transformed shared 
##D21.shared.log<- (D21.shared / rowSums(D21.shared)) * 100 # percentages

lefse.D21.shared.log<-D21.shared.log[ ,as.vector(lefse.otus$OTU)]
#filters shared file down to top 10 OTUs with highest LDA values

Igg.stat<-read.delim(file="D21.IgGposneg.no1433.1444.txt",header = F, row.names = 1)
#read in file with Igg status
lefse.D21.shared.log.meta<-merge(Igg.stat,lefse.D21.shared.log,  by= 'row.names')
row.names(lefse.D21.shared.log.meta)=lefse.D21.shared.log.meta$Row.names
lefse.D21.shared.log.meta$Row.names=NULL
lefse.D21.log.neg=lefse.D21.shared.log.meta[lefse.D21.shared.log.meta$V2=="IgG_negative",]
lefse.D21.log.neg$V2 =NULL
lefse.D21.log.pos=lefse.D21.shared.log.meta[lefse.D21.shared.log.meta$V2=="IgG_positive",]
lefse.D21.log.pos$V2 =NULL
lefse.pos<- t(lefse.D21.log.pos)
lefse.neg<- t(lefse.D21.log.neg)

#This function allows for a .taxonomy file to be converted so that it shows 
#the phylum and the last level with OTU  
make.tax<-function(taxonomy){
  new.taxonomy=taxonomy
  new.taxonomy$Phyla=NULL
  new.taxonomy$Classification_lvl100=NULL
  
  for (i in 1:length(new.taxonomy$OTU)){
    current.taxlist =  unlist(strsplit(as.character(new.taxonomy$Taxonomy[i]),');',fixed=TRUE))
    current.phyla = unlist(strsplit(as.character(current.taxlist[2]),'(',fixed=TRUE))[1]
    best="Unclassifed"
    for (j in 3:length(current.taxlist)){ 
      current.tax =  unlist(strsplit(as.character(current.taxlist[j]),'(',fixed=TRUE)) 
      if (as.numeric(current.tax[2])!=100){
        break 
      } 
      else{ 
        best = as.character(current.tax[1])
      }
    }
    current.otu=unlist(as.integer(sub("Otu","",new.taxonomy$OTU[i])))[1]
    current.otu=paste("OTU ", as.character(current.otu))
    current.phyla=gsub("_", " ", current.phyla)
    best=gsub("_", " ", best)
    new.taxonomy$Phyla[i]=current.phyla
    new.taxonomy$Classification_lvl100[i]= paste(best, " (", current.otu,")",sep ="")
  }
  #new.taxonomy$OTU=NULL
  #new.taxonomy$Size=NULL
  return(new.taxonomy)
}

tax<-read.table(file='CDIclear.final.0.03.cons.taxonomy.copy', header=TRUE)
taxa<-make.tax(taxonomy=tax)
taxa$Size<-NULL
taxa$Taxonomy<-NULL

lefse.pos.tax= merge(lefse.pos, taxa, by.x="row.names", by.y="OTU", all.x =T)
lefse.pos.tax.lda<-merge(lefse.pos.tax, lefse.otus, by.x="Row.names", by.y="OTU")
lefse.pos.tax.lda = lefse.pos.tax.lda[ order(lefse.pos.tax.lda$LDA),]
row.names(lefse.pos.tax.lda) =lefse.pos.tax.lda$Classification_lvl100
lefse.pos.tax.lda$Classification_lvl100=NULL
lefse.pos.tax.lda$Row.names=NULL
lefse.pos.tax.lda$Phyla=NULL
lefse.pos.tax.lda$LogMaxMean=NULL
lefse.pos.tax.lda$Class=NULL
lefse.pos.tax.lda$LDA=NULL
lefse.pos.tax.lda$pValue=NULL
lefse.pos.tax.lda=t(lefse.pos.tax.lda)

lefse.neg.tax= merge(lefse.neg, taxa, by.x="row.names", by.y="OTU", all.x =T)
lefse.neg.tax.lda<-merge(lefse.neg.tax, lefse.otus, by.x="Row.names", by.y="OTU")
lefse.neg.tax.lda = lefse.neg.tax.lda[ order(lefse.neg.tax.lda$LDA),]
row.names(lefse.neg.tax.lda) =lefse.neg.tax.lda$Classification_lvl100
lefse.neg.tax.lda$Classification_lvl100=NULL
lefse.neg.tax.lda$Row.names=NULL
lefse.neg.tax.lda$Phyla=NULL
lefse.neg.tax.lda$LogMaxMean=NULL
lefse.neg.tax.lda$Class=NULL
lefse.neg.tax.lda$LDA=NULL
pval.all=round(lefse.neg.tax.lda$pValue,4)
lefse.neg.tax.lda$pValue=NULL
lefse.neg.tax.lda=t(lefse.neg.tax.lda)

for (index in 1:length(pval.all)) {
  if (pval.all[index] <= 0.001) {
    pval.all[index] <- paste('= ', as.character(pval.all[index]), ' ***', sep='')
  }
  else if (pval.all[index] <= 0.01) {
    pval.all[index] <- paste('= ', as.character(pval.all[index]), ' **', sep='')
  }
  else if (pval.all[index] <= 0.05) {
    pval.all[index] <- paste('= ', as.character(pval.all[index]), ' *', sep='')
  }
  else {
    pval.all[index] <- paste('= ', as.character(pval.all[index]), ' n.s.', sep='')
  }
}

# Subset to data to cages
lefse.neg.tax.lda.c145 <- lefse.neg.tax.lda[c("1451D21","1452D21","1453D21","1454D21"),]
lefse.neg.tax.lda.c150 <- lefse.neg.tax.lda[c("1501D21","1502D21"),]
lefse.pos.tax.lda.c143 <- lefse.pos.tax.lda[c("1431D21","1432D21","1434D21"),]
lefse.pos.tax.lda.c144 <- lefse.pos.tax.lda[c("1441D21","1442D21","1443D21"),]
lefse.pos.tax.lda.c146 <- lefse.pos.tax.lda[c("1461D21","1462D21"),]
lefse.pos.tax.lda.c147 <- lefse.pos.tax.lda[c("1471D21","1472D21"),]
lefse.pos.tax.lda.c150 <- lefse.pos.tax.lda[c("1503D21","1504D21","1505D21"),]

# Reformat genera names to italics
lefsa_names <- as.list(colnames(lefse.neg.tax.lda))
lefsa_names[[10]] <- bquote(paste(italic(.('Akkermansia')), .('(OTU  3)'), sep=' '))

#Plotting
#plotting relative abundaces on log scale 
pdf(file='~/Desktop/repos/AdaptiveImmunity_and_Clearance/figures/figure_4C.pdf', width=6, height=6, useDingbats=FALSE)

par(mar=c(3,17,1,1), xaxs='r', mgp=c(2,1,0))
plot(1, type='n', ylim=c(0.8, (ncol(lefse.neg.tax.lda)*2)-0.8), xlim=c(0,4), 
     ylab='', xlab='Number of Reads (per 10000)', xaxt='n', yaxt='n', cex.lab=0.9)
box(which = "plot", lty = "solid", col ="grey80", lwd=5)
index <- 1
for(i in colnames(lefse.neg.tax.lda)){
  # Negative cages
  stripchart(at=index+0.35, lefse.neg.tax.lda.c145[,i], 
             pch=18, col='#3797a5', method='jitter', jitter=0.15, cex=1.4, lwd=1.2, add=TRUE)
  stripchart(at=index+0.35, lefse.neg.tax.lda.c150[,i], 
             pch=8, col='#3797a5', method='jitter', jitter=0.15, cex=1.2, lwd=1.2, add=TRUE)
  # Positive cages
  stripchart(at=index-0.35, lefse.pos.tax.lda.c143[,i], 
             pch=9, col='#a55637', method='jitter', jitter=0.15, cex=1.2, lwd=1.2, add=TRUE)
  stripchart(at=index-0.35, lefse.pos.tax.lda.c144[,i], 
             pch=13, col='#a55637', method='jitter', jitter=0.15, cex=1.2, lwd=1.2, add=TRUE)
  stripchart(at=index-0.35, lefse.pos.tax.lda.c146[,i], 
             pch=17, col='#a55637', method='jitter', jitter=0.15, cex=1.2, lwd=1.2, add=TRUE)
  stripchart(at=index-0.35, lefse.pos.tax.lda.c147[,i], 
             pch=15, col='#a55637', method='jitter', jitter=0.15, cex=1.2, lwd=1.2, add=TRUE)
  stripchart(at=index-0.35, lefse.pos.tax.lda.c150[,i], 
             pch=19, col='#a55637', bg='#a55637', method='jitter', jitter=0.15, cex=1, lwd=1.2, add=TRUE)
  if (i != colnames(lefse.neg.tax.lda)[length(colnames(lefse.neg.tax.lda))]){
    abline(h=index+1, lty=2)
  }
  segments(median(lefse.neg.tax.lda[,i]), index+0.6, median(lefse.neg.tax.lda[,i]), index+0.1, lwd=2.5) #adds line for median
  segments(median(lefse.pos.tax.lda[,i]), index-0.6, median(lefse.pos.tax.lda[,i]), index-0.1, lwd=2.5)
  index <- index + 2
}

axis(1, at=c(0:4), labels=c(0,10,100,1000,10000), las=1, cex.axis=0.8, col="grey80", col.ticks = "grey60")
minors <- c(0.1,0.28,0.44,0.58,0.7,0.8,0.88,0.94,0.98) 
axis(side=1, at=minors, label=rep('',length(minors)), tck=-0.01, col="grey50", col.ticks ="grey50") 
axis(side=1, at=minors+1, label=rep('',length(minors)), tck=-0.01, col="grey50", col.ticks ="grey50") 
axis(side=1, at=minors+2, label=rep('',length(minors)), tck=-0.01, col="grey50", col.ticks ="grey50") 
axis(side=1, at=minors+3, label=rep('',length(minors)), tck=-0.01,col="grey50", col.ticks ="grey50") 

legend('bottomright', legend=c('Vehicle', 'IgG Positive'),
       pch=21, pt.bg=c('#3797a5','#a55637'), bg='white', pt.cex=1.4, cex=0.75)

axis(2, at=seq(1,index-2,2)+0.5, labels=do.call(expression, lefsa_names), las=1, line=-0.5, tick=F, cex.axis=1,col="grey80", col.ticks = "grey60")
italic_p <- lapply(1:length(pval.all), function(x) bquote(paste(italic('p'), .(pval.all[x]), sep=' ')))
axis(2, at=seq(1,index-2,2)-0.5, labels=do.call(expression, italic_p), las=1, line=-0.5, tick=F, font=3,cex.axis=0.8) 

dev.off()



