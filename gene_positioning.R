############################
# Author: Anna Mathioudaki #
############################

# <><><> Required Libraries <><><>
library(argparse)

# <><><> Create Parser Argument <><><>
parser <- ArgumentParser()

parser$add_argument("-genes", "--genes", required=TRUE, help="Chromosome-Start-geneID-proteinID")
parser$add_argument("-chrom_sizes", "--chromosome_sizes", required=TRUE, help="Chromosome Sizes")
parser$add_argument("-output", "--output", required=TRUE, help="Output PDF plot")

args <- parser$parse_args()

genes <- read.table(args$genes, sep='\t', header = F) 
genes[,1] <- paste('chr', genes[,1], sep='')

b <- read.table(args$chromosome_sizes, stringsAsFactors = F, sep = "\t", header = F)
a <- as.numeric(nrow(b) -1)

#This will be my max chromosome length value
#I will use it for plotting
#I want the second largest since the largest is on Un chromosome
second_largest <- max(b[,2][b[,2]!=max(b[,2])])


pdf(args$output, height=2.7, width=9)
layout(matrix(1:nrow(b), ncol=nrow(b), byrow=F))
for(i in 1:a) {
  print(i)
  e <- b[i,1]
  print(b[i, 2])
  plot(c(0,0), c(0,b[i,2]), col='bisque4', xlab="", ylab="", ylim= c(0, second_largest), axes=F, pch=19, cex=0.7, type="l", lty=1, lwd=7, main=e)
  axis(side=2,pos=0, labels=F, tick=F)
  k <- genes[genes[,1]==e,]
  points(rep(0, length(k[,2])), k[,2], pch=4, cex=0.5)
}
dev.off()
