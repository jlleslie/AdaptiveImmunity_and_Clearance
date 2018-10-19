###Supplmental Figures for Clearance Paper 

#setworking directory
setwd("~/Desktop/AdaptiveImmunity_and_Clearance/data")

library(ggplot2)
library(grid)
library(scales)
library(vegan)
library(gtable)
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

############# Supp Figure 1A: Colonization in WT donor mice. 
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
S1A = donor.plot+ 
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
S1A = S1A + labs(y = "CFU per Gram Feces")
S1A = S1A + geom_hline(aes(yintercept=100), colour = "gray10", size = 0.9, linetype=2)
S1A = S1A + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))
S1A

############# Supp Figure 1B: Anti-C difficile toxin A titer in donor mice 
#read in the data 
antitoxin<-read.delim(file="AntitoxinA_IgGtiter_5ugmlcoat_Apri142017.txt")
# Any sample that did not have a detectable titer was given a value of 0 
# One mouse did not have sample to test and therefore got a value of NA 
#Pull out anti-toxin A IgG data for recipient mice and donor mice
recipient<-antitoxin[antitoxin$Genotype=="RAG", ]
donors<-antitoxin[antitoxin$Genotype=="WT", ]
# this can be done  using  geneotype to do this

#Statistics 
#For the donor mice the LOD for this assay was a titer of 1200
#For the purpose of stats set 0 values (ie values that no titer was detected) to LOD of 1200 divided by sqare root of 2 
fill.in.lod<-1200/sqrt(2) #fill.in.lod= 848.5281
#Replace 0 with fill in LOD 
donors$AntitoxinA_IgG_Titer<-replace(donors$AntitoxinA_IgG_Titer,donors$AntitoxinA_IgG_Titer==0,fill.in.lod)
#Testing if values from each group are from distint distributions
#Pull out Anti-toxin IgG  measured in recpient mice that got splenocytes from infected donors 
donor.infect<-donors[donors$Treatment_1=="630",]
donor.infect_AntiAigg<-c(donor.infect$AntitoxinA_IgG_Titer)
#Pull out total IgG measured in mice that got splenocytes from uninfected donors 
donor.mock<-donors[donors$Treatment_1=="mock",]
donor.mock_AntiAigg<-c(donor.mock$AntitoxinA_IgG_Titer)
wilcox.test(donor.infect_AntiAigg, donor.mock_AntiAigg, exact=F)
#data:  donor.infect_AntiAigg and donor.mock_AntiAigg
#W = 20, p-value = 0.009277
#alternative hypothesis: true location shift is not equal to 0
#No correction is required because there is only one comparison 

#to make it clear that the values that were not detected vs the values that were detected at the LOD I will replace the fill.in.lod= 848.5281 with -1000
donors$AntitoxinA_IgG_Titer<-replace(donors$AntitoxinA_IgG_Titer, donors$AntitoxinA_IgG_Titer==fill.in.lod, -10000) 

#plotting
#order the dataset so that it will plot vehicle first on the left
donors$Treatment_Grp<-factor(donors$Treatment_Grp, levels = c("uninfected_donor", "630_infected_donor"))
#assign colors
colors.don<-c("uninfected_donor"="grey", "630_infected_donor"="black")
#plot
donor.antitoxin.plot<-ggplot(donors, aes(x=Treatment_Grp, y=AntitoxinA_IgG_Titer, fill= factor(Treatment_Grp), color=factor(Treatment_Grp)))+ 
  geom_dotplot(binaxis = "y", stackdir="center", dotsize = 1.3) +
  scale_color_manual(values = rep("black",2)) +
  scale_fill_manual(values = colors.don, limits = c("uninfected_donor", "630_infected_donor")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="grey50") +
  scale_y_continuous( limits = c(-10000, 300500), labels = scales::comma, breaks = c(300000, 200000, 100000)) +
  geom_hline(aes(yintercept=1200), colour = "gray50", size = 1, linetype=2)+
  ylab("Serum Anti-TcdA IgG Titer")
S1B = donor.antitoxin.plot + 
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
    ,axis.text.y=element_text(size=11)
    ,axis.title.y=element_text(size=11)
    ,axis.title.x=element_blank()
    ,axis.text.x=element_blank()
    ,plot.margin = unit(c(1,1,2,1), "lines")
  )
S1B
#labeling the plot 
#Create text to lable plots 
gtext.doninfect<-textGrob("Infected",gp = gpar(fontsize = 10))  
gtext.donmock<-textGrob("Unifected", gp = gpar(fontsize = 10))
gtext.1200<-textGrob("1,200", gp=gpar(fontsize =11))
gtext.lod<-textGrob("(LOD)", gp=gpar(fontsize =11))
gtext.2star<-textGrob("**", gp = gpar(fontsize = 20)) 

S1B.A= S1B + annotation_custom(gtext.2star, xmin = 1.5, xmax = 1.5,  ymin = 300300, ymax = 300350) +  #adding 2 stars for comparsion between infected vs mock 
  annotate("segment", x=1, xend=2, y = 300200, yend = 300200, colour = "black", size = 0.7) +
  annotation_custom(gtext.doninfect, xmin = 2, xmax = 2, ymin = -50000, ymax = -30000) +
  annotation_custom(gtext.donmock, xmin = 1, xmax = 1, ymin  = -50000, ymax = -30000) +
  annotation_custom(gtext.1200, xmin =0.1, xmax= 0.1, ymin = 1300, ymax=1400) +
  annotation_custom(gtext.lod, xmin =0.3, xmax= 0.3, ymin = 1300, ymax=1400)

g1 = ggplotGrob(S1B.A)
g1$layout$clip[g1$layout$name=="panel"] <- "off"
grid.draw(g1)




#Supp figure 2: Colonization of other mice included in Random Forest 

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

#Supp Figure 3: Random Forest Model Using Whole Pre-treatment Community 
#please see file NEW_Figure5.R for code for this figure 




#Supp Figure 4: Relative Abundace of OTU 3(Akkermansia)

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





