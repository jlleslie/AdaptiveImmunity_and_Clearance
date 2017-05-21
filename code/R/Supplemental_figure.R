###Supplmental Figures for Clearance Paper 

#setworking directory
setwd("~/Desktop/AdaptiveImmunity_and_Clearance/data")

library(ggplot2)
library(grid)
library(scales)


# Figure S1: Colonization in WT donor mice. 
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

donor.col<-read.delim(file= "Adoptivetransfer_donor_colonization.txt", header = T)
#Replace the 100 in the mice that actually had undectable levles with LOD/squarroot of 2
#the LOD is 100 
fill.in.lod<-100/sqrt(2) #fil.in.lod = 70.71068
donor.col$CFU_g<-replace(donor.col$CFU_g,donor.col$CFU_g== 0,fill.in.lod)

##Determine the Median and IQR for CFU grouped by Cage
donor.col.med<-summaryMED(donor.col, measurevar="CFU_g", metadata=c("Cage","Day"), na.rm=TRUE)

donor.plot<-ggplot(donor.col.med, aes(x=Day, y= CFU_g, group=Cage, color=factor(Cage)))+  
  geom_point(size=6)+
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

#Supplmental figure: Colonization of other mice included in 

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
