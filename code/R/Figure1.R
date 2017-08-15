#Analysis for "The role of adaptive immunity in Clostridium difficile clearance in a murine model of infection"
# Figure 1: Results from adoptive transfer experiment

# Start with clean environment
rm(list=ls())
gc()

#load packages used for all figure 1 plots 
library(ggplot2)
library(grid)
library(gtable)

#Figure 1 A is a timeline drawn in illustrator 

setwd("~/Desktop/AdaptiveImmunity_and_Clearance/data")
#change this to where your file is located 

## Immune response to C. difficile (anti-toxin A IgG titers)
# Data used to gernate figure 1B and 1D 

#read in the data 
antitoxin<-read.delim(file="AntitoxinA_IgGtiter_5ugmlcoat_Apri142017.txt")
# Any sample that did not have a detectable titer was given a value of 0 
# One mouse did not have sample to test and therefore got a value of NA 

#Pull out anti-toxin A IgG data for recipient mice and donor mice
recipient<-antitoxin[antitoxin$Genotype=="RAG", ]
donors<-antitoxin[antitoxin$Genotype=="WT", ]
# this can be done  using  geneotype to do this


### Donor Mouse Data  Figure 2B

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
one.B = donor.antitoxin.plot + 
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
one.B

#labeling the plot 
#Create text to lable plots 
gtext.doninfect<-textGrob("Infected",gp = gpar(fontsize = 10))  
gtext.donmock<-textGrob("Unifected", gp = gpar(fontsize = 10))
gtext.1200<-textGrob("1,200", gp=gpar(fontsize =11))
gtext.lod<-textGrob("(LOD)", gp=gpar(fontsize =11))
gtext.2star<-textGrob("**", gp = gpar(fontsize = 20)) 

one.B= one.B + annotation_custom(gtext.2star, xmin = 1.5, xmax = 1.5,  ymin = 300300, ymax = 300350) +  #adding 2 stars for comparsion between infected vs mock 
  annotate("segment", x=1, xend=2, y = 300200, yend = 300200, colour = "black", size = 0.7) +
  annotation_custom(gtext.doninfect, xmin = 2, xmax = 2, ymin = -50000, ymax = -30000) +
  annotation_custom(gtext.donmock, xmin = 1, xmax = 1, ymin  = -50000, ymax = -30000) +
  annotation_custom(gtext.1200, xmin =0.1, xmax= 0.1, ymin = 1300, ymax=1400) +
  annotation_custom(gtext.lod, xmin =0.3, xmax= 0.3, ymin = 1300, ymax=1400)

g1 = ggplotGrob(one.B)
g1$layout$clip[g1$layout$name=="panel"] <- "off"
grid.draw(g1)

### Recipient Mice Data Anti-Toxin A IgG Titers 
##data for figure 1D 

#Statistics 
#For the recipient mice the LOD for this assay was a titer of 50
#For the purpose of stats set 0 values (ie values that no titer was detected) to LOD of 50 divided by sqare root of 2 
fill.in.lod<-50/sqrt(2) #fill.in.lod= 35.35534
#Replace 0 with fill in LOD 
recipient$AntitoxinA_IgG_Titer<-replace(recipient$AntitoxinA_IgG_Titer,recipient$AntitoxinA_IgG_Titer==0,fill.in.lod)
#Testing if values from each group are from distint distributions
#Pull out Anti-toxin IgG  measured in recpient mice that got splenocytes from infected donors 
rec_splen.infect<-recipient[recipient$Treatment_2=="infected_splenocytes",]
rec_splen.infect_AntiAigg<-c(rec_splen.infect$AntitoxinA_IgG_Titer)
#Pull out total IgG measured in mice that got splenocytes from uninfected donors 
rec_splen.uninfect<-recipient[recipient$Treatment_2=="mock_splenocytes",]
rec_splen.uninfect_AntiAigg<-c(rec_splen.uninfect$AntitoxinA_IgG_Titer)
#Pull out total IgG measuredIgG in mice that got vehicle 
rec_vehicle<-recipient[recipient$Treatment_2=="vehicle",]
rec_vehicl_AntiAigg<-c(rec_vehicle$AntitoxinA_IgG_Titer)

wilcox.test(rec_splen.infect_AntiAigg, rec_vehicl_AntiAigg, exact=F)
#data:  rec_splen.infect_AntiAigg and rec_vehicl_AntiAigg
#W = 48, p-value = 0.008087
#alternative hypothesis: true location shift is not equal to 0
wilcox.test(rec_splen.uninfect_AntiAigg, rec_vehicl_AntiAigg, exact=F)
#data:  rec_splen.uninfect_AntiAigg and rec_vehicl_AntiAigg
#W = 18, p-value = NA
#alternative hypothesis: true location shift is not equal to 0
wilcox.test(rec_splen.infect_AntiAigg, rec_splen.uninfect_AntiAigg,exact=F)
#data:  rec_splen.infect_AntiAigg and rec_splen.uninfect_AntiAigg
#W = 48, p-value = 0.008087
#alternative hypothesis: true location shift is not equal to 0

#Correcting P-values for mutiple comparisons
recip_pvals<-c(0.008087,0.008087, NA)
round(p.adjust(recip_pvals, method = "BH"),3)
# Antitoxin A titer in mice  that recived splenocytes from infected donors vs mice that got vehicle: 0.008
# Antitoxin A titer in mice that recived splenocytes from infected donors vs uninfected donors:  0.008 


#Plotting the data 
#to make it clear that the values that were not detected vs the values that were detected at the LOD I will replace the fill.in.lod= 35.35534 with -500
recipient$AntitoxinA_IgG_Titer<-replace(recipient$AntitoxinA_IgG_Titer, recipient$AntitoxinA_IgG_Titer==fill.in.lod, -500) 
#order the dataset so that it will plot vehicle first on the left
recipient$Treatment_2<-factor(recipient$Treatment_2, levels = c("vehicle", "mock_splenocytes", "infected_splenocytes"))

#assign colors to different treatment groups
colors.recip<-c("infected_splenocytes"="#f91780", "mock_splenocytes"= "#fa8c17", "vehicle"="#0095a3")
#plot 
recipient.antitoxin.plot <-ggplot(recipient, aes(x=Treatment_2, y=AntitoxinA_IgG_Titer,  fill= factor(Treatment_2), colour= factor(Treatment_2)))+ 
  geom_dotplot(binaxis = "y", stackdir="center", dotsize = 1.3) +
  scale_color_manual(values = rep("black",3)) +
  scale_fill_manual(values = colors.recip, limits = c("infected_splenocytes", "vehicle")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="grey50") +
  scale_y_continuous( limits = c(-500, 4500)) +
  geom_hline(aes(yintercept=50), colour = "gray50", size = 1, linetype=2)+
  ylab(" Serum Anti-Toxin A IgG Titer")
one.D = recipient.antitoxin.plot + 
  #eliminates background, gridlines and key border and other changes to theme
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
one.D #plots what you have so far

# Add in the lables for treatment groups outside of the plot
gtext.spleninfect<-textGrob("Splenocytes\n(infected donor)", gp = gpar(fontsize = 10))  
gtext.splenmock<-textGrob("Splenocytes\n(uninfected donor)", gp = gpar(fontsize = 10))  
gtext.vhe<-textGrob("Vehicle",gp = gpar(fontsize = 10)) 
gtext.star<-textGrob("*", gp = gpar(fontsize = 20)) 
gtext.ns <-textGrob("ns", gp = gpar(fontsize = 13)) 

one.D = one.D+ annotation_custom(gtext.2star, xmin=2, xmax=2, ymin=4225, ymax=4275) +  #adding 2 stars for comparsion between splenocytes-infect vs vehicle 
  annotate("segment", x = 1, xend = 3, y = 4225, yend = 4225, colour = "black", size = 0.7) +
  annotation_custom(gtext.2star, xmin=2.5, xmax=2.5, ymin=4500, ymax=4550) +  #adding 2 stars for comparsion between splenocytes- uninfect vs vehicle 
  annotate("segment", x = 2, xend = 3, y = 4500, yend = 4500, colour = "black", size = 0.7) +
  annotation_custom(gtext.spleninfect, xmin=3, xmax=3,  ymin=-1200, ymax=-1000) +
  annotation_custom(gtext.splenmock, xmin=2, xmax=2,  ymin=-1200, ymax=-1000) +
  annotation_custom(gtext.vhe, xmin=1, xmax=1,ymin=-1200, ymax=-1000)

g2 = ggplotGrob(one.D)
g2$layout$clip[g2$layout$name=="panel"] <- "off"
grid.draw(g2)




#Figure 1D
## Measurement of Total IgG in mice at time of harvest 
#read in the data 
IgG<-read.delim(file="TotalIgG_April_5_2017.txt", header = T)

#Statistics 
#Testing if values from each group are from distint distributions
#For the purpose of stats set 0 values (ie values that no IgG was detected) to LOD of 1.56 divided by sqare root of 2 
fill.in.lod<-1.56/sqrt(2) #fill.in.lod= 1.103087
#Replace 0 with fill in LOD 
IgG$Total_IgG<-replace(IgG$Total_IgG, IgG$Total_IgG==0,fill.in.lod)
#Pull out total IgG measured in mice that got splenocytes from infected donors 
spl_infected_totIgg<-IgG[IgG$Treatment_2=="infected_splenocytes",]
spl_infected_totIgg<-c(spl_infected_totIgg$Total_IgG)
#Pull out total IgG measured in mice that got splenocytes from uninfected donors 
spl_uninfected_totIgg<-IgG[IgG$Treatment_2=="mock_splenocytes",]
spl_uninfected_totIgg<-c(spl_uninfected_totIgg$Total_IgG)
#Pull out total IgG measuredIgG in mice that got vehicle 
vehicle_totIgg<-IgG[IgG$Treatment_2=="vehicle",]
vehicle_totIgg<-c(vehicle_totIgg$Total_IgG)

wilcox.test(spl_infected_totIgg, vehicle_totIgg, exact=F)
#data:  spl_infected_totIgg and vehicle_totIgg
#W = 51, p-value = 0.003508
#alternative hypothesis: true location shift is not equal to 0
wilcox.test(spl_uninfected_totIgg, vehicle_totIgg, exact=F)
#data:  spl_uninfected_totIgg and vehicle_totIgg
#W = 33, p-value = 0.009622
#alternative hypothesis: true location shift is not equal to 0
wilcox.test(spl_infected_totIgg, spl_uninfected_totIgg,exact=F)
#W = 24.5, p-value = 0.8135

#Correcting P-values for mutiple comparisons
totigg_pvals<-c(0.003508,0.009622,0.8135)
round(p.adjust(totigg_pvals, method = "BH"),3)
# Total IgG in infected splenocytes vs vehicle: 0.011
# Total IgG in uninfected splenocytes ve vehicle :  0.014
# Total IgG in infected splenocytes vs uninfected splenocytes groups: 0.814

#Plotting the data 
#to make it clear that the values that were not 


#Plot data 

#to clearly show points that were not detected (vs detected at LOD), change values to be below LOD line in the figure 
#fill.in.lod = 1.032376
#Because I had previously replaced 0 with the fill.in.lod I have to change that value to -500
IgG$Total_IgG<-replace(IgG$Total_IgG, IgG$Total_IgG==fill.in.lod, -500) 
#if you want to change the values again you will have to replace fill.in.lod with -500 etc. 

#assign colors to different treatment groups
colors<-c("infected_splenocytes"="#f91780", "mock_splenocytes"= "#fa8c17", "vehicle"="#0095a3")

##order the dataset so that it will plot vehicle first on the left
IgG$Treatment_2<-factor(IgG$Treatment_2, levels = c("vehicle", "mock_splenocytes", "infected_splenocytes"))

#Make a jitter plot of the data 
igg.plot<-ggplot(IgG, aes(x=Treatment_2, y=Total_IgG, fill=factor(Treatment_2)))+ 
  geom_jitter(aes(shape=21) , width = 0.2, height = 0.01, size = 4)+
  scale_shape_identity()+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="grey50") +
  geom_hline(aes(yintercept=1.56), colour = "grey50", size = 1, linetype=2) +
  scale_fill_manual(values = colors) +
  scale_y_continuous( limits = c(-1000, 9000)) 
 
one.C= igg.plot + 
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
one.C  = one.C + labs(y = "Serum IgG ng/ml")
one.C 

# Add in the lables for treatment groups outside of the plot
gtext.spleninfect<-textGrob("Splenocytes\n(infected donor)", gp = gpar(fontsize = 10))  
gtext.splenmock<-textGrob("Splenocytes\n(uninfected donor)", gp = gpar(fontsize = 10))  
gtext.vhe<-textGrob("Vehicle",gp = gpar(fontsize = 10)) 
gtext.star<-textGrob("*", gp = gpar(fontsize = 20)) 
gtext.ns <-textGrob("ns", gp = gpar(fontsize = 13)) 


one.C  = one.C  + annotation_custom(gtext.star, xmin=2, xmax=2, ymin=8150, ymax=8200) +  #adding single star for comparsion between splenocytes-infect vs vehicle 
          annotate("segment", x = 1, xend = 3, y = 8130, yend = 8130, colour = "black", size = 0.7) +
          annotation_custom(gtext.star, xmin=1.5, xmax=1.5, ymin=8700, ymax=8750) +  #adding single star for comparsion between splenocytes- uninfect vs vehicle 
          annotate("segment", x = 1, xend = 2, y = 8625, yend = 8650, colour = "black", size = 0.7) +
          annotation_custom(gtext.ns, xmin=2.5, xmax=2.5, ymin=9200, ymax=9300) + #adding ns for comparsion between splenocytes- uninfect vs plenocytes- uninfect
          annotate("segment", x = 2, xend =3, y = 8925, yend = 8950, colour = "black", size = 0.7)+
           annotation_custom(gtext.spleninfect, xmin=3, xmax=3,  ymin=-2400, ymax=-2000) +
          annotation_custom(gtext.splenmock, xmin=2, xmax=2, ymin=-2400, ymax=-2000) +
          annotation_custom(gtext.vhe, xmin=1, xmax=1, ymin=-2400, ymax=-2000)

g3 = ggplotGrob(one.C)
g3$layout$clip[g$layout$name=="panel"] <- "off"
grid.draw(g3)




