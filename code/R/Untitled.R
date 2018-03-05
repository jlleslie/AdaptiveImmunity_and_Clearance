
#Figure 4D  CFU for at end of experiment by genotype and group
#Colonization overtime
cfutime_data<-read.table(file='Colonization_Overtime_630_Allexperiments_copy.txt', header=TRUE)
#Note the data is reported such that if no colonies were seen on the 10^-2 plate, that sample was reported as having 100 CFU/g feces ie the LOD 

##2014 Exeperiment data, RAG and WT mice
cfutime_data<-cfutime_data[cfutime_data$Experiment == "2014",]

#Pulls out Mice that were infected 
cfutime_data.infected<-cfutime_data[cfutime_data$Treatment_1 =="630",]
cfutime_data.infectedD40<-cfutime_data.infected[cfutime_data.infected$Day=="40",]

fill.in.lod<-100/sqrt(2) #fil.in.lod = 70.71068
cfutime_data.infectedD40$CFU_g<-replace(cfutime_data.infectedD40$CFU_g,cfutime_data.infectedD40$CFU_g==100,fill.in.lod)
col=c("16"="#61bfee", "978"="#61bfee","18"="#bb5fa1","977" ="#bb5fa1")
shape_A=as.numeric(c("RAG"="1", "WT" = "19"))

D40.plot<-ggplot(cfutime_data.infectedD40, aes(x= factor(Cage), y= CFU_g, group=Cage, color=factor(Cage), shape=Genotype))+  
  geom_jitter(size=6, width = 0.1)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="grey50") +
  scale_shape_manual(values=shape_A) +
  scale_color_manual(values=col)


#theme with white background
d= D40.plot + 
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
d1 = d+ labs(y = " CFU per Gram Feces")
d2 = d1+ geom_hline(aes(yintercept=100), colour = "gray10", size = 0.9, linetype=2)
d3 = d2 +  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))
d3