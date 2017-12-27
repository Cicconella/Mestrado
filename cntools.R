source("http://bioconductor.org/biocLite.R")
library("CNTools")
library(stargazer)

samples = LETTERS[1:6]

a = c("A",5,10,25,25-10,5)
a = rbind(a,c("A",5,40,55,15, 2))
b = c("B",5,15,35,20,5)
c = c("C",5,5,45,40,5)
e = c("E",5,5,20,15,6)
f = c("F",5,45,55,10,1)


aux = as.data.frame(rbind(a,b,c,e,f))
aux
colnames(aux) = c("ID", "chrom", "loc.start", "loc.end", "num.mark",  "seg.mean")
aux$ID = as.character(aux$ID)
aux$chrom = as.numeric(as.character(aux$chrom))
aux$loc.start = as.numeric(as.character(aux$loc.start))
aux$loc.end = as.numeric(as.character(aux$loc.end))
aux$num.mark = as.numeric(as.character(aux$num.mark))
aux$seg.mean = as.numeric(as.character(aux$seg.mean))

aux

stargazer(aux, rownames = F, summary = F)



attach(aux)

n=length(samples)
aux <- data.frame(c(aux$ID,samples), c(aux$chrom,rep(unique(aux$chrom), n)), 
                  c(aux$loc.start,rep(1, n)), c(aux$loc.end,rep(max(aux$loc.end), n)), 
                  c(aux$num.mark,rep(1, n)), c(aux$seg.mean,rep(-1, n)))
colnames(aux) = c("ID", "chrom", "loc.start", "loc.end", "num.mark",  "seg.mean")
aux

aux$ID = as.character(aux$ID)
str(aux)
head(aux)
dim(aux)

seg <-  CNSeg(aux)
seg
rsByregion <- getRS(seg, by = "region", imput = TRUE, XY = FALSE, what = "max")
cnvA = rs(rsByregion)
cnvA

dim(cnvA)
cnvA = cbind(cnvA[,1:3],apply(cnvA[,-c(1:3)],2,arruma))
cnvA

stargazer(cnvA, rownames = F, summary = F)



