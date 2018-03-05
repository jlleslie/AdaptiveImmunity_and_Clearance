#Analysis for "The gut microbiota is associated with clearance of Clostridium difficile infection independent of adaptive immunity"
# Figure 2: Results from adoptive transfer experiment

# Start with clean environment
rm(list=ls())
gc()

#load packages used for all figure 1 plots 
deps <- c('ggplot2', 'vegan', 'reshape2', 'grid', 'scales', 'gtable');
for (dep in deps){
  if (dep %in% installed.packages()[,'Package'] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

setwd("~/Desktop/AdaptiveImmunity_and_Clearance/data")
#change this to where your file is located 


####### Figure 2A
##Total Serum IgG in recpeient mice  
IgG<-read.delim(file="TotalIgG_April_5_2017.txt", header = T)
 
#Statistics 
#The LOD for this assay was a titer of 1.56
#For the purpose of stats set 0 values (ie values that no titer was detected) to LOD of 50 divided by sqare root of 2 
fill.in.lod<-1.56/sqrt(2) 
#fill.in.lod= 1.103
IgG$Total_IgG<-replace(IgG$Total_IgG, IgG$Total_IgG==0,1.103)
#Replace 0 with fill in LOD 
PT = pairwise.wilcox.test(IgG$Total_Ig, 
                          IgG$Treatment_2, 
                          p.adjust.method="BH")

PT
#data:  IgG$Total_Ig and IgG$Treatment_2 
                  #infected_splenocytes mock_splenocytes
#mock_splenocytes        0.814                -               
#  vehicle               0.011                0.014     

#pairwise wilxcox tests with Benjamini & Hochberg correction 
#Results 
  #Infected vs mock: NS
  # Infected vs vehicle: 0.011
  # Vehicle  vs mock : 0.014
#plotting the data 
#to clearly show points that were not detected (vs detected at LOD), change values to be below LOD line in the figure 
IgG$Total_IgG<-replace(IgG$Total_IgG, IgG$Total_IgG==1.103, -500) 
#assign colors to different treatment groups
colors<-c("infected_splenocytes"="#f91780", "mock_splenocytes"= "#fa8c17", "vehicle"="#0095a3")
##order the dataset so that it will plot vehicle first on the left
IgG$Treatment_2<-factor(IgG$Treatment_2, levels = c("vehicle", "mock_splenocytes", "infected_splenocytes"))
shape_cage<-c("143"= 21, "144"=22, "145" =23, "146"=24, "147" =13, "150" =9 , "150A" =8)
#find median of each group
infected_splenocytes.igg.med<- median(IgG[IgG$Treatment_2=="infected_splenocytes",11]) 
mock_splenocytes.igg.med<- median(IgG[IgG$Treatment_2=="mock_splenocytes",11]) 
vehicle.igg.med<- median(IgG[IgG$Treatment_2=="vehicle",11]) 

#Make a jitter plot of the data 
igg.plot<-ggplot(IgG, aes(x=Treatment_2, y=Total_IgG, color=factor(Treatment_2), fill=factor(Treatment_2),shape=factor(Cage)))+ 
  scale_shape_manual(values= shape_cage)+
  geom_jitter(width = 0.35, height = 0.01, size= 5)+
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(-1000, 9000), labels = scales::comma, breaks = c(8000, 6000, 4000, 2000)) +
  geom_hline(aes(yintercept=1.56), colour = "grey50", size = 1, linetype=2)
fig2A= igg.plot + 
 theme (
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
fig2A  = fig2A+ labs(y = "Serum IgG ng/ml")
fig2A
# Add in the lables for treatment groups outside of the plot
gtext.spleninfect<-textGrob("Splenocytes\n(infected donor)", gp = gpar(fontsize = 10))  
gtext.splenmock<-textGrob("Splenocytes\n(uninfected donor)", gp = gpar(fontsize = 10))  
gtext.vhe<-textGrob("Vehicle",gp = gpar(fontsize = 10)) 
gtext.star<-textGrob("*", gp = gpar(fontsize = 20)) 
gtext.ns <-textGrob("ns", gp = gpar(fontsize = 13)) 
gtext.1<-textGrob("1.56", gp=gpar(fontsize =11))
gtext.lod<-textGrob("(LOD)", gp=gpar(fontsize =11))

fig2A = fig2A + annotation_custom(gtext.star, xmin=2, xmax=2, ymin=8150, ymax=8200) +  #adding single star for comparsion between splenocytes-infect vs vehicle 
  annotate("segment", x = 1, xend = 3, y = 8130, yend = 8130, colour = "black", size = 0.7) +
  annotation_custom(gtext.star, xmin=1.5, xmax=1.5, ymin=8700, ymax=8750) +  #adding single star for comparsion between splenocytes- uninfect vs vehicle 
  annotate("segment", x = 1, xend = 2, y = 8625, yend = 8650, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=2.5, xmax=2.5, ymin=9200, ymax=9300) + #adding ns for comparsion between splenocytes- uninfect vs splenocytes- infect
  annotate("segment", x = 2, xend =3, y = 8925, yend = 8950, colour = "black", size = 0.7)+
  annotate("segment", x=0.8, xend=1.2, y=vehicle.igg.med, yend=vehicle.igg.med, color="grey50", size=1) +
  annotate("segment", x=1.8, xend=2.2, y=infected_splenocytes.igg.med, yend=infected_splenocytes.igg.med, color="grey50", size=1) +
  annotate("segment", x=2.8, xend=3.2, y=mock_splenocytes.igg.med, yend=mock_splenocytes.igg.med, color="grey50", size=1)+ 
  #the three segments above draw in median lines for each group 
   annotation_custom(gtext.spleninfect, xmin=3, xmax=3,  ymin=-2400, ymax=-2000) +
  annotation_custom(gtext.splenmock, xmin=2, xmax=2, ymin=-2400, ymax=-2000) +
  annotation_custom(gtext.vhe, xmin=1, xmax=1, ymin=-2400, ymax=-2000)+
  annotation_custom(gtext.1, xmin =0.06, xmax= 0.08, ymin = 1.56, ymax=3) +
  annotation_custom(gtext.lod, xmin =0.25, xmax= 0.25, ymin = 1.56, ymax=3)

g3 = ggplotGrob(fig2A)
g3$layout$clip[g3$layout$name=="panel"] <- "off"
grid.draw(g3)


############Figure 2B
### Recipient Mice Data Anti-Toxin A IgG Titers 

#read in the data 
antitoxin<-read.delim(file="AntitoxinA_IgGtiter_5ugmlcoat_Apri142017.txt")
# Any sample that did not have a detectable titer was given a value of 0 
# One mouse did not have sample to test and therefore got a value of NA 

#Pull out anti-toxin A IgG data for recipient mice and donor mice
recipient<-antitoxin[antitoxin$Genotype=="RAG", ]
# this can be done  using  geneotype 

#Statistics 
#For the recipient mice the LOD for this assay was a titer of 50
#For the purpose of stats set 0 values (ie values that no titer was detected) to LOD of 50 divided by sqare root of 2 
fill.in.lod<-50/sqrt(2) 
#fill.in.lod= 35.35534
#Replace 0 with fill in LOD 
recipient$AntitoxinA_IgG_Titer<-replace(recipient$AntitoxinA_IgG_Titer,recipient$AntitoxinA_IgG_Titer==0,35.35534)
#Testing if values from each group are from distint distributions
#Pull out Anti-toxin IgG  measured in recpient mice that got splenocytes from infected donors 
PT = pairwise.wilcox.test(recipient$AntitoxinA_IgG_Titer, 
                          recipient$Treatment_2, 
                          p.adjust.method="BH")

PT
                      #infected_splenocytes mock_splenocytes
#mock_splenocytes        0.0081                 -               
#  vehicle               0.0081                 -    

#pairwise wilxcox tests with Benjamini & Hochberg correction 
#Results 
#Infected vs mock:  0.0081
# Infected vs vehicle: 0.0081
# Vehicle  vs mock : - (both are same)


#Plotting the data 
#to make it clear that the values that were not detected vs the values that were detected at the LOD I will replace the fill.in.lod= 35.35534 with -500
recipient$AntitoxinA_IgG_Titer<-replace(recipient$AntitoxinA_IgG_Titer, recipient$AntitoxinA_IgG_Titer==35.35534, -500) 
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
  scale_y_continuous( limits = c(-500, 4500), breaks = c(4000, 3000, 2000, 1000)) +
  geom_hline(aes(yintercept=50), colour = "gray50", size = 1, linetype=2)+
  ylab("Serum Anti-TcdA IgG Titer")
fig2B = recipient.antitoxin.plot + 
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
fig2B #plots what you have so far

# Add in the lables for treatment groups outside of the plot
gtext.spleninfect<-textGrob("Splenocytes\n(infected donor)", gp = gpar(fontsize = 10))  
gtext.splenmock<-textGrob("Splenocytes\n(uninfected donor)", gp = gpar(fontsize = 10))  
gtext.vhe<-textGrob("Vehicle",gp = gpar(fontsize = 10)) 
gtext.2star<-textGrob("**", gp = gpar(fontsize = 20)) 
gtext.ns <-textGrob("ns", gp = gpar(fontsize = 13)) 
gtext.50<-textGrob("50", gp=gpar(fontsize =11))
gtext.lod<-textGrob("(LOD)", gp=gpar(fontsize =11))

fig2B = fig2B + annotation_custom(gtext.2star, xmin=2, xmax=2, ymin=4225, ymax=4275) +  #adding 2 stars for comparsion between splenocytes-infect vs vehicle 
  annotate("segment", x = 1, xend = 3, y = 4225, yend = 4225, colour = "black", size = 0.7) +
  annotation_custom(gtext.2star, xmin=2.5, xmax=2.5, ymin=4500, ymax=4550) +  #adding 2 stars for comparsion between splenocytes- uninfect vs vehicle 
  annotate("segment", x = 2, xend = 3, y = 4500, yend = 4500, colour = "black", size = 0.7) +
  annotation_custom(gtext.spleninfect, xmin=3, xmax=3,  ymin=-1200, ymax=-1000) +
  annotation_custom(gtext.splenmock, xmin=2, xmax=2,  ymin=-1200, ymax=-1000) +
  annotation_custom(gtext.vhe, xmin=1, xmax=1,ymin=-1200, ymax=-1000) +
  annotation_custom(gtext.50, xmin =0.08, xmax= 0.08, ymin = 50, ymax=100) +
  annotation_custom(gtext.lod, xmin =0.25, xmax= 0.25, ymin = 50, ymax=100)

g2 = ggplotGrob(fig2B)
g2$layout$clip[g2$layout$name=="panel"] <- "off"
grid.draw(g2)


#Figure 2C
#Colonization levels over time 

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

#read in the data 
colonization<-read.table(file="Colonization_RAG_Recips.txt", header=T)
#LOD is 100 CFU/g feces
fill.in.lod<-100/sqrt(2)
colonization$CFU_g[colonization$CFU_g == "0"] <-fill.in.lod
#changes 0 values to LoD/sqt(2) so that it is clear they were not detected 

colonization.noD1<-colonization[colonization$Day != "1", ]
#removes Day 1 data from the data set, the data from this day is comprommised 
#because the lables that were placed on the tubes threw off the weights and they were diluted way too much (and don't have a mechanism to fix it)

cfu.treat<-summaryMED(colonization.noD1, measurevar="CFU_g", metadata=c("Treatment_2","Day"), na.rm=TRUE)
cfu.cage<-summaryMED(colonization.noD1, measurevar="CFU_g", metadata=c("Cage","Day"), na.rm=TRUE)
group.cols<-c("infected_splenocytes"="#f91780", "mock_splenocytes"= "#fa8c17", "vehicle"="#0095a3",
              "143"="#f91780", "146"="#f91780" , "150"="#f91780", "144"= "#db9204", "147" ="#db9204", "145"="#006670", "150A"= "#006670")
#Plot data
cfu.plot<-ggplot(cfu.treat, aes(x=Day, y=CFU_g, color=factor(Treatment_2)))+ 
  geom_line(size=2.5)+
  geom_line(data=cfu.cage, aes(x=Day, y=CFU_g, color=factor(Cage)), size= 0.55,linetype = 2) +
  scale_color_manual(values=group.cols)+
  scale_x_continuous(breaks = c(3,6,9,12,15,18,21,24), limits =c(3, 26))
#theme with white background
fig2C= cfu.plot + 
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
fig2C.1 = fig2C + labs(y = " CFU per Gram Feces", x = "Day Post Infection")
fig2C.2 = fig2C.1+ geom_hline(aes(yintercept=100), colour = "gray50", size = 1, linetype=2)
fig2C.3 =fig2C.2 + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) 
fig2C.3
# Add in the lables for treatment groups outside of the plot
gtext.spleninfect<-textGrob("Splenocytes\n(infected donor)", gp = gpar(fontsize = 8))  
gtext.splenmock<-textGrob("Splenocytes\n(uninfected donor)", gp = gpar(fontsize = 8))  
gtext.vhe<-textGrob("Vehicle",gp = gpar(fontsize = 8)) 
fig2C.3  = fig2C.3 + 
  annotate("rect", xmin = 17, xmax = 18, ymin = 100000, ymax = 200000, fill = "#0095a3")+  
  annotation_custom(gtext.vhe, xmin=20, xmax=20,  ymin = 5, ymax =5.3 ) +
  annotate("rect", xmin = 17, xmax = 18, ymin = 10000, ymax = 20000, fill = "#fa8c17")+  
  annotation_custom(gtext.splenmock, xmin = 21, xmax = 21,  ymin = 4, ymax = 4.3) +
  annotate("rect", xmin = 17, xmax = 18, ymin = 1000, ymax = 2000, fill = "#f91780")+  
  annotation_custom(gtext.spleninfect, xmin = 21, xmax = 21,  ymin =3, ymax = 3.3) 
g3 = ggplotGrob(fig2C.3)
g3$layout$clip[g3$layout$name=="panel"] <- "off"
grid.draw(g3)




###Figure 2d
## Colonication at the end of the experiment 

## Pull out CFU at time at D26 (23 days post trainsfer of splenocytes)
colonization.D26<-colonization[colonization$Day == "26", ]

#Replace the 0 in the mice that had undectable levles with LOD/squarroot of 2
#the LOD is 100 
fill.in.lod<-100/sqrt(2) #fil.in.lod = 70.71068
colonization.D26$CFU_g<-replace(colonization.D26$CFU_g,colonization.D26$CFU_g==0,fill.in.lod)

#Testing if values from each group are from distint distributions
pairwise.wilcox.test(colonization.D26$CFU_g, 
                          colonization.D26$Treatment_2, 
                          p.adjust.method="BH")

#data:  colonization.D26$CFU_g and colonization.D26$Treatment_2 

#                   infected_splenocytes mock_splenocytes
#mock_splenocytes       1.00                 -               
#vehicle              1.00                 0.72            
#P value adjustment method: BH 
# CFU at D26 Infected vs Vehicle :1
# CFU at D26 Infected vs Mock: 1
# CFU at D26 Mock vs Vehicle:  0.72

#find median of each group
infected_splenocytes.cfu.med<- median(colonization.D26[colonization.D26$Treatment_2=="infected_splenocytes",11]) 
mock_splenocytes.cfu.med<- median(colonization.D26[colonization.D26$Treatment_2=="mock_splenocytes",11]) 
vehicle.cfu.med<- median(colonization.D26[colonization.D26$Treatment_2=="vehicle",11]) 

##Plotting D26 Data 
colors<-c("infected_splenocytes"="#f91780", "mock_splenocytes"= "#fa8c17", "vehicle"="#0095a3")
#order the dataset so that it will plot vehicle first on the left
colonization.D26$Treatment_2<-factor(colonization.D26$Treatment_2, levels = c("vehicle", "mock_splenocytes", "infected_splenocytes"))
shape_cage<-c("143"= 21, "144"=22, "145" =23, "146"=24, "147" =13, "150" =9 , "150A" =8)
d26.plot<-ggplot(colonization.D26, aes(x=Treatment_2, y=CFU_g, color=factor(Treatment_2), fill=factor(Treatment_2), shape=factor(Cage)))+ 
  #stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="grey50") +
  scale_shape_manual(values= shape_cage)+
  geom_jitter(width = 0.2, height = 0.01, size= 5)+
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)), limits=c(10,4*10^8)) +
  labs(y = expression(paste(" CFU ", "per Gram Feces")))+
  geom_hline(aes(yintercept=100), colour = "gray50", size = 1, linetype=2)
#theme with white background
fig2D = d26.plot + 
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
fig2D  = fig2D  + annotation_custom(gtext.ns, xmin=2, xmax=2, ymin = 8.1, ymax= 8.2) + 
  annotate("segment", x = 1, xend = 3, y = 100000000, yend = 100000000, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=2.5, xmax=2.5, ymin=8.35, ymax=8.6) +  
  annotate("segment", x = 2, xend = 3, y = 200000000, yend = 200000000, colour = "black", size = 0.7) +
  annotation_custom(gtext.ns, xmin=1.5, xmax=1.5, ymin=8.6, ymax=8.8) +
  annotate("segment", x = 1, xend = 2, y =300000000, yend = 300000000, colour = "black", size = 0.7)+
  annotate("segment", x=0.8, xend=1.2, y=vehicle.cfu.med, yend=vehicle.cfu.med, color="grey50", size=1) +
  annotate("segment", x=1.8, xend=2.2, y=infected_splenocytes.cfu.med, yend=infected_splenocytes.cfu.med, color="grey50", size=1) +
  annotate("segment", x=2.8, xend=3.2, y=mock_splenocytes.cfu.med, yend=mock_splenocytes.cfu.med, color="grey50", size=1)+ 
  annotation_custom(gtext.spleninfect, xmin=3, xmax=3,  ymin=0, ymax=0) +
  annotation_custom(gtext.splenmock, xmin=2, xmax=2,  ymin=0, ymax=0) +
  annotation_custom(gtext.vhe, xmin=1, xmax=1, ymin=0, ymax=0)
g = ggplotGrob(fig2D)
g$layout$clip[g$layout$name=="panel"] <- "off"
grid.draw(g)




