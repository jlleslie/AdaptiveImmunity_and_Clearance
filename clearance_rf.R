
# Start with clean environment
rm(list=ls())
gc()
set.seed(8619)

# Load dependencies
deps <- c('randomForest', 'vegan', 'gtools');
for (dep in deps){
  if (dep %in% installed.packages()[,'Package'] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

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

random_Forest <- function(training_data, feat){

  attach(training_data)
  # ntree and mtry calculations from Segal et al. (2004)
  num1 <- round(length(as.character(rownames(data[which(training_data[,feat]==as.vector(levels(training_data[,feat]))[1]),]))* 0.623))
  num2 <- round(length(as.character(rownames(data[which(training_data[,feat]==as.vector(levels(training_data[,feat]))[2]),]))* 0.623))
  rf_factor <- max(c(round(num1/num2), round(num2/num2))) * 3 
  trees <- round(length(colnames(training_data)) - 1) * rf_factor
  tries <- round(sqrt(length(colnames(training_data)) - 1))
  
  # Run random forest
  rf_object <- randomForest(training_data[,feat]~., data=training_data, importance=TRUE, replace=FALSE, do.trace=500, err.rate=TRUE, ntree=trees, mtry=tries)
  
  rm(num1, num2, rf_factor, trees, tries)
  detach(training_data)
  return(rf_object)
}

# Define data files
metadata <- '~/Desktop/Repositories/clearance_2017/Adaptiveimmuneclear_metadata_noD40.42.txt'
shared <- '~/Desktop/Repositories/clearance_2017/CDIclear.final.shared'
taxonomy <- '~/Desktop/Repositories/clearance_2017/CDIclear.final.0.03.cons.taxonomy'

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

#--------------------------------------------------------------------#

# Decide on optimal subsample level
sub_size <- floor(as.numeric(quantile(rowSums(shared), probs=0.09)))

# Plot abundances to pick a subsample size
pdf(file='~/Desktop/Repositories/clearance_2017/subsample_plot.pdf', width=8, height=4)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))
par(mar=c(3,3,1,1), mgp=c(2,1,0))
plot(sort(rowSums(shared)), pch=20, xlab='Sample', ylab='Abundance')
abline(h=sub_size, lty=4, col='red', lwd=2) # Show the cutoff
text(x=200, y=85000, 'Subsample Size:', cex=0.8)
text(x=200, y=80000, sub_size)
mtext('A', side=2, line=2, las=2, adj=1, padj=-8, cex=1.5)

# Rarefaction for random sample to test size
rarecurve(shared[736,], sample=sub_size)
mtext('B', side=2, line=2, las=2, adj=1, padj=-8, cex=1.5)
dev.off()

#--------------------------------------------------------------------#

# Format data
# Rarefy and filter shared file
shared <- as.data.frame(t(shared))
shared <- shared[, colSums(shared) > sub_size] # Loses 7% samples
for (index in 1:ncol(shared)){
  shared[,index] <- t(rrarefy(shared[,index], sample=sub_size))}
rm(index)
shared <- as.data.frame(t(shared))
shared <- filter_table(shared) # Loses 5007 OTUs

# Subset taxonomy to remaining OTUs
taxonomy <- subset (taxonomy, rownames(taxonomy) %in% colnames(shared))

# Merge datasets
shared <- merge(metadata, shared, by='row.names', all.x=FALSE, all.y=FALSE) # Drops 339 samples
rownames(shared) <- shared$Row.names
shared$Row.names <- NULL
rm(metadata)

# Subset to groups of interest for analysis
# Cleared vs Colonized
cleared_colonized <- subset(shared, Colonization630 != 'uncolonized')
cleared_colonized <- subset(cleared_colonized, Day %in% c(-15,-12,1)) # needs work...

cleared_colonized$Co_Housed <- NULL
cleared_colonized$Day <- NULL
cleared_colonized$Treatment_2 <- NULL

# Adoptive transfer
adoptive_transfer <- subset(shared, Treatment_2 %in% c('infected_splenocytes','mock_splenocytes'))
# further subsetting?
adoptive_transfer$Day <- NULL
adoptive_transfer$Co_Housed <- NULL
adoptive_transfer$Colonization630 <- NULL

# What outcome based on co-housing is of interest?
cohousing <- subset(shared, Colonization630 != 'uncolonized')
cohousing$Colonization630 <- NULL
cohousing$Co_Housed <- NULL
cohousing$Treatment_Grp <- NULL
cohousing$Treatment_1 <- NULL
cohousing$Treatment_2 <- NULL

rm(shared)

#--------------------------------------------------------------------#

# Get perumations of groups to equalize random forest
options(expressions=1e5)
test <- combinations(n=length(rownames(cleared_colonized[which(cleared_colonized$Colonization630==as.vector(levels(cleared_colonized$Colonization630))[2]),])), 
                     r=length(rownames(cleared_colonized[which(cleared_colonized$Colonization630==as.vector(levels(cleared_colonized$Colonization630))[1]),])), 
                     v=rownames(cleared_colonized[which(cleared_colonized$Colonization630==as.vector(levels(cleared_colonized$Colonization630))[2]),]))
# n -> size of source vector
# r -> size of target vector
# v -> source vector


# Run random forest iteratively for all permutations





RF.importances <- importance(data.randomForest, type=1)
sig.sort.final <- subset(RF.importances, RF.importances > abs(min(RF.importances)))





#--------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/rf_plot.pdf'
pdf(file=plot_file, width=6, height=12)
layout(matrix(c(1,2,
                3,4,
                5,6), nrow=3, ncol=2, byrow=TRUE))

#--------------------------------------------------------------------#

# Cleared vs colonized
# RF mean decrease accuracy
par(mar=c(3,3,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(importances$Metabolite_score, labels=importances$Compound_name,
         lcolor=NA, cex=1.2, groups=importances$abx, color='black',
         xlab='Importance Score', xlim=c(0,10), pch=19, lwd=3,
         gcolor=c('darkmagenta',wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5],'forestgreen'))
mtext('A', side=2, line=2, las=2, adj=1, padj=-17, cex=1.7)
segments(x0=rep(0, 15), y0=c(1:5, 8:11, 14:18, 21:25, 28:32), 
         x1=rep(12, 15), y1=c(1:5, 8:11, 14:18, 21:25, 28:32), lty=2) # Dotted lines

# OTU abundance changes
plot(1, type='n', ylim=c(0,length(clinda_otus)*2), xlim=c(0,4), 
     ylab='', xlab=expression(paste('Normalized Abundance (',log[10],')')), xaxt='n', yaxt='n')
title('Clindamycin-pretreated', line=0.5, cex.main=1.1, font.main=2, col.main=clinda_col)
index <- 1
for(i in colnames(cleared_otu)){
  stripchart(at=index-0.35, jitter(cleared_otu[,i], amount=1e-5), 
             pch=21, bg='firebrick1', method='jitter', jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  stripchart(at=index+0.35, jitter(colonized_otu[,i], amount=1e-5), 
             pch=21, bg='dodgerblue1', method='jitter', jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  if (i != colnames(clinda_mock_otu)[length(colnames(clinda_mock_otu))]){
    abline(h=index+1, lty=2)
  }
  segments(median(clinda_mock_otu[,i]), index-0.4, median(clinda_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(clinda_infected_otu[,i]), index+0.4, median(clinda_infected_otu[,i]), index, lwd=2)
  index <- index + 2
}
axis(1, at=seq(0,4,1), label=c('0','10','100','1000','10000'))
minor.ticks.axis(1, 10, mn=0, mx=4)
legend('topright', legend=c('Cleared', 'Colonized'),
       pch=c(21, 21), pt.bg=c('firebrick1','dodgerblue1'), bg='white', pt.cex=1.4, cex=0.9)
formatted <- lapply(1:length(clinda_otus), function(i) bquote(paste(italic(.(clinda_genera[i])), .(clinda_otus[i]), sep=' ')))
axis(2, at=seq(1,index-2,2)+0.4, labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
axis(2, at=seq(1,index-2,2), labels=clinda_phyla, las=1, line=-0.5, tick=F, cex.axis=1.1) 
formatted <- lapply(1:length(clinda_pvalues), function(i) bquote(paste(italic('p'), .(clinda_pvalues[i]), sep=' ')))
axis(2, at=seq(1,index-2,2)-0.5, labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 

#---------------------#

# Adoptive transfer

#---------------------#

# Cohousing


dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

