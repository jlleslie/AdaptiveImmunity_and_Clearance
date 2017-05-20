
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

# Filter out columns that have values in at least 3 samples (ignores first column if needed)
filter_table <- function(data) {
  drop <- c()
  if (class(data[,1]) != 'character') {
    if (sum(data[,1] != 0) < 3) {
      drop <- c(drop, colnames(data)[1])
    }
  }
  for (index in 2:ncol(data)) {
    if (sum(data[,index] != 0) < 3) {
      drop <- c(drop, colnames(data)[index])
    }
  }
  filtered_data <- data[,!(colnames(data) %in% drop)]
  return(filtered_data)
}

# Define data files
metadata <- '~/Desktop/Repositories/clearance_2017/data/Adaptiveimmuneclear_metadata_noD40.42.txt'
shared <- '~/Desktop/Repositories/clearance_2017/data/Adaptiveimmuneclear_noD40.42.0.03.filter.0.03.subsample.shared'
taxonomy <- '~/Desktop/Repositories/clearance_2017/data/clearance.formatted.taxonomy'

# Define output file
plot_file <- '~/Desktop/Repositories/clearance_2017/figures/day0_RF_plot.pdf'

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
# Decide on optimal subsample level
#sub_size <- floor(as.numeric(quantile(rowSums(shared), probs=0.09)))

# Rarefy and filter shared file
#shared <- as.data.frame(t(shared))
#shared <- shared[, colSums(shared) > sub_size] # Loses 7% samples
#for (index in 1:ncol(shared)){
#  shared[,index] <- t(rrarefy(shared[,index], sample=sub_size))}
#rm(index, sub_size)
#shared <- as.data.frame(t(shared))
#shared <- filter_table(shared) # Loses 5007 OTUs

# Subset taxonomy to remaining OTUs
taxonomy <- subset(taxonomy, rownames(taxonomy) %in% colnames(shared))

# Merge datasets
shared <- merge(metadata, shared, by='row.names', all.x=FALSE, all.y=FALSE) # Drops 339 samples
rownames(shared) <- shared$Row.names
shared$Row.names <- NULL
rm(metadata)

# Subset to groups of interest for analysis
# Cleared vs Colonized
cleared_colonized <- subset(shared, Colonization630 != 'uncolonized')
cleared_colonized$Colonization630 <- factor(cleared_colonized$Colonization630)
cleared_colonized_day0 <- subset(cleared_colonized, Day == 0)
cleared_colonized_day0$Co_Housed <- NULL
cleared_colonized_day0$Day <- NULL
cleared_colonized_day0$Treatment_2 <- NULL
rm(cleared_colonized)

rm(shared)

#--------------------------------------------------------------------#

# Determine optimal ntree and mtry
factor1 <- as.vector(levels(cleared_colonized_day0$Colonization630))[1]
factor2 <- as.vector(levels(cleared_colonized_day0$Colonization630))[2]
num1 <- round(length(which(cleared_colonized_day0$Colonization630 == factor1)) * 0.623)
num2 <- round(length(which(cleared_colonized_day0$Colonization630 == factor2)) * 0.623)
ntree_multiplier <- max(c((num1/num2), (num2/num2))) * 3 
trees <- round(ncol(cleared_colonized_day0) - 1) * ntree_multiplier
tries <- round(sqrt(ncol(cleared_colonized_day0) - 1))
rm(factor1, factor2, num1, num2, ntree_multiplier)

# Run random forest and assess predictive value
cleared_day0_rf <- randomForest(cleared_colonized_day0$Colonization630~., 
                                data=cleared_colonized_day0, importance=TRUE, replace=FALSE, 
                                do.trace=FALSE, err.rate=TRUE, ntree=trees, mtry=tries)
rm(trees, tries)

# Retreive importance and overall error rate
day0_importances <- importance(cleared_day0_rf, type=1)
day0_accuracy <- paste('Accuracy = ',as.character((100-round((median(cleared_day0_rf$err.rate[,1])*100), 2))),'%',sep='')

# Merge remaining important OTUs with taxonomy
day0_importances <- merge(day0_importances, taxonomy, by='row.names', all.x=TRUE)
rownames(day0_importances) <- day0_importances$Row.names
day0_importances$Row.names <- NULL
colnames(day0_importances)[1] <- 'MDA'
rm(taxonomy)

# Subset to the most important OTUs and sort
day0_importances <- subset(day0_importances, day0_importances$MDA > abs(min(day0_importances$MDA)))
day0_importances <- as.data.frame(day0_importances[order(-day0_importances$MDA),])
day0_importances <- day0_importances[1:10,]
day0_importances <- as.data.frame(day0_importances[order(day0_importances$MDA),])
rm(cleared_day0_rf)

# Subset important OTU abundances from shared file
cleared_day0_shared <- cleared_colonized_day0[, c(1,which(colnames(cleared_colonized_day0) %in% rownames(day0_importances)))]
rm(cleared_colonized_day0)

# Break into experimental groups and match order to importance
colonized_day0_shared <- subset(cleared_day0_shared, Colonization630 == 'colonized')
colonized_day0_shared$Colonization630 <- NULL
colonized_day0_shared <- colonized_day0_shared[,match(rownames(day0_importances), colnames(colonized_day0_shared))] 
cleared_day0_shared <- subset(cleared_day0_shared, Colonization630 == 'cleared')
cleared_day0_shared$Colonization630 <- NULL
cleared_day0_shared <- cleared_day0_shared[,match(rownames(day0_importances), colnames(cleared_day0_shared))] 

#--------------------------------------------------------------------#

# Testing for significant change abundance
day0_pvalues <- c()
for (index in 1:ncol(colonized_day0_shared)){
  day0_pvalues[index] <- wilcox.test(cleared_day0_shared[,index], colonized_day0_shared[,index], exact=FALSE)$p.value
}
day0_pvalues <- round(p.adjust(day0_pvalues, method='BH'), 3)
for (index in 1:length(day0_pvalues)) {
  if (day0_pvalues[index] == 0) {
    day0_pvalues[index] <- '< 0.001 ***'
  }
  else if (day0_pvalues[index] == 0.001) {
    day0_pvalues[index] <- paste('= ', as.character(day0_pvalues[index]), ' ***', sep='')
  }
  else if (day0_pvalues[index] <= 0.01) {
    day0_pvalues[index] <- paste('= ', as.character(day0_pvalues[index]), ' **', sep='')
  }
  else if (day0_pvalues[index] <= 0.05) {
    day0_pvalues[index] <- paste('= ', as.character(day0_pvalues[index]), ' *', sep='')
  }
  else {
    day0_pvalues[index] <- paste('= ', as.character(day0_pvalues[index]), ' n.s.', sep='')
  }
}
day0_importances$pvalues <- day0_pvalues
rm(day0_pvalues)

# Log10 transform for easier viewing
cleared_day0_shared <- log10(cleared_day0_shared + 1)
colonized_day0_shared <- log10(colonized_day0_shared + 1)

#--------------------------------------------------------------------#

# Set up plotting environment
pdf(file=plot_file, width=11, height=6)
layout(matrix(c(1,2,2), nrow=1, ncol=3, byrow=TRUE))

#---------------------#

# Day 0
# RF mean decrease accuracy
par(mar=c(1.8,3,1,1), xaxs='i', xaxt='n', xpd=FALSE, mgp=c(2,0.2,0))
dotchart(day0_importances$MDA, labels=rownames(day0_importances),
         lcolor=NA, cex=1.2, color='black', 
         xlab='', xlim=c(5,13), pch=19, lwd=3)
segments(x0=rep(5, 10), y0=c(1:10), x1=rep(13, 10), y1=c(1:10), lty=2) # Dotted lines
legend('bottomright', legend=day0_accuracy, pt.cex=0, cex=1.2, bty='n')
par(xaxt='s')
axis(side=1, at=c(5:13), labels=c(0,6:13), cex.axis=1.2, tck=-0.025)
axis.break(1, 5.5, style='slash')
mtext('Mean Decrease Accuracy', side=1, padj=1.8, cex=0.9)
mtext('A', side=2, line=2, las=2, adj=1, padj=-13.2, cex=1.7)

# OTU abundance differences
par(mar=c(3,20,1,1), xaxs='r', mgp=c(2,1,0))
plot(1, type='n', ylim=c(0.8, (ncol(cleared_day0_shared)*2)-0.8), xlim=c(0,4), 
     ylab='', xlab='Abundance', xaxt='n', yaxt='n', cex.lab=1.4)
index <- 1
for(i in colnames(cleared_day0_shared)){
  stripchart(at=index+0.35, cleared_day0_shared[,i], 
             pch=21, bg='firebrick1', method='jitter', jitter=0.15, cex=1.7, lwd=0.5, add=TRUE)
  stripchart(at=index-0.35, colonized_day0_shared[,i], 
             pch=21, bg='dodgerblue1', method='jitter', jitter=0.15, cex=1.7, lwd=0.5, add=TRUE)
  if (i != colnames(cleared_day0_shared)[length(colnames(cleared_day0_shared))]){
    abline(h=index+1, lty=2)
  }
  segments(median(cleared_day0_shared[,i]), index+0.6, median(cleared_day0_shared[,i]), index+0.1, lwd=2.5) #adds line for median
  segments(median(colonized_day0_shared[,i]), index-0.6, median(colonized_day0_shared[,i]), index-0.1, lwd=2.5)
  index <- index + 2
}
axis(side=1, at=c(0:4), label=c('0','10','100','1000','10000'), cex.axis=1.2, tck=-0.02)
minors <- c(0.1,0.28,0.44,0.58,0.7,0.8,0.88,0.94,0.98)
axis(side=1, at=minors, label=rep('',length(minors)), tck=-0.01)
axis(side=1, at=minors+1, label=rep('',length(minors)), tck=-0.01)
axis(side=1, at=minors+2, label=rep('',length(minors)), tck=-0.01)
axis(side=1, at=minors+3, label=rep('',length(minors)), tck=-0.01)
legend('topright', legend=c('Cleared', 'Colonized'),
       pch=c(21, 21), pt.bg=c('firebrick1','dodgerblue1'), bg='white', pt.cex=1.7, cex=1.2)
axis(2, at=seq(1,index-2,2)+0.6, labels=rownames(day0_importances), las=1, line=-0.5, tick=F, cex.axis=1.4)
formatted_taxa <- lapply(1:nrow(day0_importances), function(x) bquote(paste(.(day0_importances$phylum[x]),'; ',italic(.(day0_importances$genus[x])), sep='')))
axis(2, at=seq(1,index-2,2), labels=do.call(expression, formatted_taxa), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
italic_p <- lapply(1:length(day0_importances$pvalues), function(x) bquote(paste(italic('p'), .(day0_importances$pvalues[x]), sep=' ')))
axis(2, at=seq(1,index-2,2)-0.6, labels=do.call(expression, italic_p), las=1, line=-0.5, tick=F, cex.axis=1.2, font=3) 
mtext('B', side=2, line=2, las=2, adj=13, padj=-13, cex=1.7)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

