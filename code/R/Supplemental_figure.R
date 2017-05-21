###Supplmental Figures for Clearance Paper 

#setworking directory
setwd("~/Desktop/AdaptiveImmunity_and_Clearance/data")
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
S1 = donor.plot+ 
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
S1 = S1+ labs(y = "CFU per Gram Feces")
S2 = S1+ geom_hline(aes(yintercept=100), colour = "gray10", size = 0.9, linetype=2)
S3 = S2 + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))
S3
