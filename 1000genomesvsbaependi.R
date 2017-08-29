## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
#biocLite("quantsmooth")
require("quantsmooth")

chr1Genomes = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/1000genomes/c1", 
                         header = T)

head(chr1Genomes)


# prepareGenomePlot example
library(quantsmooth)
# construct genomic positions
CHR<-sample(1,1,replace=TRUE)  # Chromosomes
MapInfo<-lengthChromosome(CHR,"bases")*runif(length(CHR)) # position on chromosome
chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE, organism="hsa")
# Chrompos returns a matrix with the positions of the elements on the plot
# You can use all kinds of base graphics functions to annotate the chromosomes
points(chrompos[,2],chrompos[,1]+0.1,pch="x",col="red")
# Show connection between 3rd and 4th element
segments(chrompos[3,2],chrompos[3,1],chrompos[4,2],chrompos[4,1],col="blue",lwd=2)

j=1
cnv = read.table(paste("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/cromo", j, sep=""), header = T)

cnv=cnv[,c(1:3)]
head(cnv)
head(chr1Genomes)
