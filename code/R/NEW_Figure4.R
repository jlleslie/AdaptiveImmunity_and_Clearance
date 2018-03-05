###Analysis for Figure 1: 
#Co-housing of RAG and WT mice 
#Question can microbiota alone enable clearance 
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

setwd("~/Desktop/AdaptiveImmunity_and_Clearance/data")


##Figure 4 Finding OTUS that discriminate IgG+ vs IgG- (vehicle) mice
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
D21.shared.log<- log10(D21.shared + 1)
#calcs Log10 transformed shared 

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


#Plotting
#plotting relative abundaces on log scale 
par(mar=c(3,17,1,1), xaxs='r', mgp=c(2,1,0))
plot(1, type='n', ylim=c(0.8, (ncol(lefse.neg.tax.lda)*2)-0.8), xlim=c(0,4), 
     ylab='', xlab='Relative Abundance', xaxt='n', yaxt='n', cex.lab=1.2)
box(which = "plot", lty = "solid", col ="grey80", lwd=5)
index <- 1
for(i in colnames(lefse.neg.tax.lda)){
  stripchart(at=index+0.35, lefse.neg.tax.lda[,i], 
             pch=21, bg='grey', method='jitter', jitter=0.15, cex=1.4, lwd=0.5, add=TRUE)
  stripchart(at=index-0.35, lefse.pos.tax.lda[,i], 
             pch=21, bg='grey20', method='jitter', jitter=0.15, cex=1.4, lwd=0.5, add=TRUE)
  if (i != colnames(lefse.neg.tax.lda)[length(colnames(lefse.neg.tax.lda))]){
    abline(h=index+1, lty=2)
  }
  segments(median(lefse.neg.tax.lda[,i]), index+0.6, median(lefse.neg.tax.lda[,i]), index+0.1, lwd=2.5) #adds line for median
  segments(median(lefse.pos.tax.lda[,i]), index-0.6, median(lefse.pos.tax.lda[,i]), index-0.1, lwd=2.5)
  index <- index + 2
}

axis(side=1, at=c(0:4), label=c('0','0.01','1','10',"100"), cex.axis=1.2, tck=-0.02, col="grey50", col.ticks ="grey50")
minors <- c(0.1,0.28,0.44,0.58,0.7,0.8,0.88,0.94,0.98)
axis(side=1, at=minors, label=rep('',length(minors)), tck=-0.01, col="grey50", col.ticks ="grey50")
axis(side=1, at=minors+1, label=rep('',length(minors)), tck=-0.01, col="grey50", col.ticks ="grey50")
axis(side=1, at=minors+2, label=rep('',length(minors)), tck=-0.01, col="grey50", col.ticks ="grey50")
axis(side=1, at=minors+3, label=rep('',length(minors)), tck=-0.01,col="grey50", col.ticks ="grey50")

legend('bottomright', legend=c('Vehicle', 'IgG Positive'),
       pch=c(21, 21), pt.bg=c('grey','grey20'), bg='white', pt.cex=1.4, cex=0.5)
axis(2, at=seq(1,index-2,2)+0.5, labels=colnames(lefse.neg.tax.lda), las=1, line=-0.5, tick=F, cex.axis=1,col="grey80", col.ticks = "grey60")
italic_p <- lapply(1:length(pval.all), function(x) bquote(paste(italic('p'), .(pval.all[x]), sep=' ')))
axis(2, at=seq(1,index-2,2)-0.5, labels=do.call(expression, italic_p), las=1, line=-0.5, tick=F, font=3,cex.axis=0.8) 





