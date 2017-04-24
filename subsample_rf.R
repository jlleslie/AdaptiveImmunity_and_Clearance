
# Start with clean environment
rm(list=ls())
gc()
set.seed(8619)

# Load dependencies
deps <- c('vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,'Package'] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep)

# Finds index of median in a vector
which.median = function(x) {
  if (length(x) %% 2 != 0) {
    which(x == median(x))
  } else if (length(x) %% 2 == 0) {
    a = sort(x)[c(length(x)/2, length(x)/2+1)]
    c(which(x == a[1]), which(x == a[2]))
  }
}

# Define data files
shared <- '~/Desktop/Repositories/clearance_2017/CDIclear.final.shared'

# Read in data and eliminate extra columns
shared <- read.delim(shared, sep='\t', header=T, row.names=2)
shared$numOtus <- NULL
shared$label <- NULL

#--------------------------------------------------------------------#

# Decide on optimal subsample level
sub_size <- floor(as.numeric(quantile(rowSums(shared), probs=0.09)))

# Find sample with median diversity
median_diversity <- names(which.median(diversity(shared, index='invsimpson', MARGIN=1)))

# Plot abundances to pick a subsample size
pdf(file='~/Desktop/Repositories/clearance_2017/subsample_plot.pdf', width=8, height=4)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))
par(mar=c(3,3,1,1), mgp=c(2,1,0))
plot(sort(rowSums(shared)), pch=20, xlab='Sample', ylab='Abundance')
abline(h=sub_size, lty=4, col='red', lwd=2) # Show the cutoff
text(x=200, y=85000, 'Subsample Size:', cex=0.8)
text(x=200, y=80000, sub_size)
mtext('A', side=2, line=2, las=2, adj=1, padj=-8, cex=1.5)

# Rarefaction for sample with median inv. Simpson diversity
rarecurve(shared[median_diversity,], sample=sub_size)
mtext('B', side=2, line=2, las=2, adj=1, padj=-8, cex=1.5)
dev.off()

#--------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

