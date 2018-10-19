
# Start with clean environment
rm(list=ls())
gc()
set.seed(8619)

# Load dependencies
deps <- c('randomForest', 'vegan', 'plotrix');
for (dep in deps){
  if (dep %in% installed.packages()[,'Package'] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep)

# Define data files
#setwd("~/Desktop/AdaptiveImmunity_and_Clearance/data")
setwd("~/Desktop/repos/clearance_2017/data")
metadata <- 'Adaptiveimmuneclear_metadata_noD40.42.tsv'
shared <- 'Adaptiveimmuneclear_noD40.42.0.03.subsample.0.03.filter.shared'
taxonomy <- 'clearance.formatted.taxonomy'


# Read in data and eliminate extra columns
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metadata$Mouse <- NULL
metadata$Cage <- NULL
metadata$Gender <- NULL
metadata$Genotype <- NULL
metadata$Year <- NULL
metadata$Treatment_Grp <- NULL
metadata$Treatment_1 <- NULL
shared <- read.delim(shared, sep='\t', header=T, row.names=2)
shared$numOtus <- NULL
shared$label <- NULL
taxonomy <- read.delim(taxonomy, sep='\t', header=T, row.names=1)
taxonomy$Size <- NULL
taxonomy$phylum <- gsub('_', ' ', taxonomy$phylum)
taxonomy$family <- gsub('_', ' ', taxonomy$family)
taxonomy$genus <- gsub('_', ' ', taxonomy$genus)

#--------------------------------------------------------------------#

# Format data

# Subset taxonomy to remaining OTUs
taxonomy <- subset(taxonomy, rownames(taxonomy) %in% colnames(shared))

# Merge datasets
shared <- merge(metadata, shared, by='row.names', all.x=FALSE, all.y=FALSE) # Drops 339 samples
rownames(shared) <- shared$Row.names
shared$Row.names <- NULL
rm(metadata)

# Subset to groups of interest for analysis
# Cleared vs Colonized
cleared_colonized <- subset(shared, shared$Colonization_stat630 != 'unchallenged')
cleared_colonized$Colonization_stat630 <- factor(cleared_colonized$Colonization_stat630)
cleared_colonized_preabx <- subset(cleared_colonized, Day %in% c(-15,-12))
cleared_colonized_preabx$Co_Housed <- NULL
cleared_colonized_preabx$Day <- NULL
cleared_colonized_preabx$Treatment_2 <- NULL
rm(cleared_colonized)
rm(shared)

# Subset to most informative OTUs from previous RF analysis
#pruned_shared <- cleared_colonized_preabx[, c('Colonization_stat630','Otu0052', 'Otu0093', 'Otu0026')]
pruned_shared <- cleared_colonized_preabx[, c('Colonization_stat630','Otu0052', 'Otu0093')]

#--------------------------------------------------------------------#

# preabx time point
# Determine optimal ntree and mtry
factor1 <- as.vector(levels(cleared_colonized_preabx$Colonization_stat630))[1]
factor2 <- as.vector(levels(cleared_colonized_preabx$Colonization_stat630))[2]
num1 <- round(length(which(cleared_colonized_preabx$Colonization_stat630 == factor1)) * 0.623)
num2 <- round(length(which(cleared_colonized_preabx$Colonization_stat630 == factor2)) * 0.623)
ntree_multiplier <- max(c((num1/num2), (num2/num2))) * 3 
trees <- round(ncol(cleared_colonized_preabx) - 1) * ntree_multiplier
tries <- round(sqrt(ncol(cleared_colonized_preabx) - 1))
rm(factor1, factor2, num1, num2, ntree_multiplier)

# Run random forest and assess predictive value
temp <- cleared_colonized_preabx
test <- temp$Colonization_stat630
temp$Colonization_stat630 <- NULL
cleared_preabx_rf <- randomForest(test~., data=temp, importance=TRUE, replace=FALSE, 
                                  do.trace=FALSE, err.rate=TRUE, ntree=trees, mtry=tries)
rm(temp, test, trees, tries)
preabx_accuracy <- paste('Accuracy = ',as.character( (100-round((tail(cleared_preabx_rf$err.rate[,1], n=1)*100), 2))),'%',sep='')
#preabx_accuracy <- 'Accuracy = 76.92%'
print(cleared_preabx_rf)
print(preabx_accuracy)

# Retreive importance and overall error rate
preabx_importances <- importance(cleared_preabx_rf, type=1)

# Get accuracy of binning from pruned feature set
factor1 <- as.vector(levels(pruned_shared$Colonization_stat630))[1]
factor2 <- as.vector(levels(pruned_shared$Colonization_stat630))[2]
num1 <- round(length(which(pruned_shared$Colonization_stat630 == factor1)) * 0.623)
num2 <- round(length(which(pruned_shared$Colonization_stat630 == factor2)) * 0.623)
ntree_multiplier <- max(c((num1/num2), (num2/num2))) * 3 
trees <- round(ncol(pruned_shared) - 1) * ntree_multiplier
tries <- round(sqrt(ncol(pruned_shared) - 1))
rm(factor1, factor2, num1, num2, ntree_multiplier)

temp <- pruned_shared
test <- as.factor(temp$Colonization_stat630)
temp$Colonization_stat630 <- NULL
pruned_rf <- randomForest(test~., data=temp, importance=TRUE, replace=FALSE, 
                                  do.trace=FALSE, err.rate=TRUE, ntree=trees, mtry=tries)
rm(temp, test, trees, tries)
pruned_accuracy <- paste('Accuracy = ',as.character( (100-round((tail(pruned_rf$err.rate[,1], n=1)*100), 2))),'%',sep='')
print(pruned_rf)
print(pruned_accuracy)

# Merge remaining important OTUs with taxonomy
preabx_importances <- merge(preabx_importances, taxonomy, by='row.names', all.x=TRUE)
rownames(preabx_importances) <- preabx_importances$Row.names
preabx_importances$Row.names <- NULL
colnames(preabx_importances)[1] <- 'MDA'

# Subset to the most important OTUs and sort
preabx_importances <- subset(preabx_importances, preabx_importances$MDA > abs(min(preabx_importances$MDA)))
preabx_importances <- as.data.frame(preabx_importances[order(-preabx_importances$MDA),])
preabx_importances <- preabx_importances[1:10,]
preabx_importances <- as.data.frame(preabx_importances[order(preabx_importances$MDA),])
rm(cleared_preabx_rf)

# Subset important OTU abundances from shared file
cleared_preabx_shared <- cleared_colonized_preabx[, c(1,which(colnames(cleared_colonized_preabx) %in% rownames(preabx_importances)))]
rm(cleared_colonized_preabx)

# Break into experimental groups
colonized_preabx_shared <- subset(cleared_preabx_shared, Colonization_stat630 == 'colonized')
colonized_preabx_shared$Colonization_stat630 <- NULL
colonized_preabx_shared <- colonized_preabx_shared[,match(rownames(preabx_importances), colnames(colonized_preabx_shared))] 
cleared_preabx_shared <- subset(cleared_preabx_shared, Colonization_stat630 == 'cleared')
cleared_preabx_shared$Colonization_stat630 <- NULL
cleared_preabx_shared <- cleared_preabx_shared[,match(rownames(preabx_importances), colnames(cleared_preabx_shared))]

#--------------------------------------------------------------------#

# Testing for significant change abundance
preabx_pvalues <- c()
for (index in 1:ncol(colonized_preabx_shared)){
  preabx_pvalues[index] <- wilcox.test(cleared_preabx_shared[,index], colonized_preabx_shared[,index], exact=FALSE)$p.value
}
preabx_pvalues <- round(p.adjust(preabx_pvalues, method='holm'), 3)
for (index in 1:length(preabx_pvalues)) {
  if (preabx_pvalues[index] == 0) {
    preabx_pvalues[index] <- '< 0.001 ***'
  }
  else if (preabx_pvalues[index] == 0.001) {
    preabx_pvalues[index] <- paste('= ', as.character(preabx_pvalues[index]), ' ***', sep='')
  }
  else if (preabx_pvalues[index] <= 0.01) {
    preabx_pvalues[index] <- paste('= ', as.character(preabx_pvalues[index]), ' **', sep='')
  }
  else if (preabx_pvalues[index] <= 0.05) {
    preabx_pvalues[index] <- paste('= ', as.character(preabx_pvalues[index]), ' *', sep='')
  }
  else {
    preabx_pvalues[index] <- paste('= ', as.character(preabx_pvalues[index]), ' n.s.', sep='')
  }
}
preabx_importances$pvalues <- preabx_pvalues
rm(preabx_pvalues)

# Log10 transform for easier viewing
cleared_preabx_shared <- log10(cleared_preabx_shared + 1)
colonized_preabx_shared <- log10(colonized_preabx_shared + 1)

#--------------------------------------------------------------------#

tiff(filename='~/Desktop/test.tiff', width=600, height= 400)


# Set up plotting environment
layout(matrix(c(1,2,2), nrow=1, ncol=3, byrow=TRUE))

#---------------------#

# Cleared vs colonized
# RF median decrease accuracy
par(mar=c(2,3,1,1), xaxs='i', xaxt='n', xpd=FALSE, mgp=c(2,0.2,0))
dotchart(preabx_importances$MDA, labels=toupper(rownames(preabx_importances)),
         lcolor=NA, cex=1.1, color='black', 
         xlab='', xlim=c(2,14), pch=19, lwd=3)
segments(x0=rep(2, 10), y0=c(1:10), x1=rep(14, 10), y1=c(1:10), lty=2) # Dotted lines
legend('bottomright', legend=preabx_accuracy, pt.cex=0, cex=1.2, bty='n')
par(xaxt='s')
axis(side=1, at=c(2,4,6,8,10,12,14), labels=c(0,4,6,8,10,12,14), cex.axis=1.2, tck=-0.025)
axis.break(1, 3, style='slash')
mtext('Mean Decrease Accuracy', side=1, padj=1.8, cex=0.9)
mtext('A', side=2, line=2, las=2, adj=1, padj=-17.3, cex=1.7)

# OTU abundance differences
# Abundance (per 10000 sequences)
par(mar=c(3,20,1,1), xaxs='r', mgp=c(2,1,0))
plot(1, type='n', ylim=c(0.8, (ncol(cleared_preabx_shared)*2)-0.8), xlim=c(0,4), 
     ylab='', xlab='Relative Abundance %', xaxt='n', yaxt='n', cex.lab=1.4)
index <- 1
for(i in colnames(cleared_preabx_shared)){
  stripchart(at=index+0.35, cleared_preabx_shared[,i], 
             pch=21, bg='deeppink', method='jitter', jitter=0.15, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=index-0.35, colonized_preabx_shared[,i], 
             pch=21, bg='darkblue', method='jitter', jitter=0.15, cex=2, lwd=0.5, add=TRUE)
  if (i != colnames(cleared_preabx_shared)[length(colnames(cleared_preabx_shared))]){
    abline(h=index+1, lty=2)
  }
  segments(median(cleared_preabx_shared[,i]), index+0.6, median(cleared_preabx_shared[,i]), index+0.1, lwd=2.5) #adds line for median
  segments(median(colonized_preabx_shared[,i]), index-0.6, median(colonized_preabx_shared[,i]), index-0.1, lwd=2.5)
  index <- index + 2
}
axis(side=1, at=c(0:4), label=c('0','0.1','1','10','100'), cex.axis=1.2, tck=-0.02)
minors <- c(0.1,0.28,0.44,0.58,0.7,0.8,0.88,0.94,0.98)
axis(side=1, at=minors, label=rep('',length(minors)), tck=-0.01)
axis(side=1, at=minors+1, label=rep('',length(minors)), tck=-0.01)
axis(side=1, at=minors+2, label=rep('',length(minors)), tck=-0.01)
axis(side=1, at=minors+3, label=rep('',length(minors)), tck=-0.01)
legend('topright', legend=c('Cleared', 'Colonized'),
      pch=c(21, 21), pt.bg=c('deeppink','darkblue'), bg='white', pt.cex=2, cex=1.1)
axis(2, at=seq(1,index-2,2)+0.6, labels=toupper(rownames(preabx_importances)), las=1, line=-0.5, tick=F, cex.axis=1.4)
formatted_taxa <- lapply(1:nrow(preabx_importances), function(x) bquote(paste(.(preabx_importances$phylum[x]),'; ',italic(.(preabx_importances$genus[x])), sep='')))
axis(2, at=seq(1,index-2,2), labels=do.call(expression, formatted_taxa), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
italic_p <- lapply(1:length(preabx_importances$pvalues), function(x) bquote(paste(italic('p'), .(preabx_importances$pvalues[x]), sep=' ')))
axis(2, at=seq(1,index-2,2)-0.6, labels=do.call(expression, italic_p), las=1, line=-0.5, tick=F, cex.axis=1.2, font=3) 
mtext('B', side=2, line=2, las=2, adj=13, padj=-17, cex=1.7)

dev.off()



#-------------------------------------------------------------------------------------------------------------------------------------#



# Plots abundances of top 2 OTUS

abund <- c(log10(pruned_shared$Otu0052 + 1), log10(pruned_shared$Otu0093 + 1))
colz <- c(pruned_shared$Colonization_stat630, pruned_shared$Colonization_stat630)
# cleared = 1
# colonized = 2
otu <- c(rep(1, nrow(pruned_shared)),rep(2, nrow(pruned_shared)))
# Otu0052 = 1
# Otu0093 = 2
pruned <- cbind(otu, colz, abund)




pdf(file='~/Desktop/repos/AdaptiveImmunity_and_Clearance/figures/Final_OTUs.pdf', width=6, height=4)
par(mar=c(3,3,1,1), xaxs='i', xaxt='n', xpd=FALSE, mgp=c(2,0.75,0), las=1)
boxplot(abund~colz*otu, data=pruned, at=c(1.25,2.05,3.25,4.05),
        col=c('white','gray'), outline=FALSE, ylim=c(0,3),
        boxlwd=2, staplewex=0, whisklwd=2, lty=1, ylab='Relative Abundance (Log10)')
legend('topleft', legend=c('Cleared','Colonized'), pt.bg=c('white','gray'), 
       pch=22, cex=1.2, pt.cex=2, box.col='white')
segments(x0=c(1.25,3.25), y0=c(2,2.8),
      x1=c(2.05,4.05), y1=c(2,2.8),
      lwd=1.5)
text(x=c(1.65,3.65), y=c(2.15,2.95), c('*','*'), font=2, cex=2)
mtext(c('OTU 52','OTU 93'), side=1, at=c(1.6,3.6), padj=0.5, cex=1.4)
mtext(c('(MDA: 0.097)','(MDA: 0.175)'), side=1, at=c(1.6,3.6), padj=2.5, cex=0.9)
box(lwd=3, col='gray')
dev.off()


#-------------------------------------------------------------------------------------------------------------------------------------#


