for(i in 1:22){
  
  cnvA = cnv[Chr==i,]
  dim(cnvA)
  head(cnvA)
  aux = cnvA[order(cnvA$Start),]
  dim(aux)
  head(aux)
  
  aux[1:20,]
  
  summary(aux[,8])
  which(aux[,8]==0)
  
  bla = aux[c(1,4,100,489),]
  
  aux = bla
  
  # Change matrix to data.frame and add an auxiliar region from the start of the 
  # chromossome
  
  
  aux <- data.frame(c(aux$Sample,aux$Sample), c(aux$Chr,rep(i, length(aux$Sample))), 
                    c(aux$Start,rep(1, length(aux$Sample))), c(aux$End,rep(max(aux$End), length(aux$Sample))), 
                    c(aux$Number,rep(1, length(aux$Sample))), c(aux$State,rep(-1, length(aux$Sample))))
  
  colnames(aux) = c("ID", "chrom", "loc.start", "loc.end", "num.mark",  "seg.mean")
  str(aux)
    aux
  dim(aux)
  
  seg <-  CNSeg(aux)
  seg
  rsByregion <- getRS(seg, by = "region", imput = TRUE, XY = FALSE, what = "max")
  cnvA = rs(rsByregion)
  head(cnvA)
  dim(cnvA)
  cnvA
  
  cnvA = apply(cnvA,2,arruma)

  cnvA
    
  chromosomes[[i]] = cnvA
}

arruma = function(x){
  x[which(x==1)] = 0
  x[which(x==2)] = 1
  x[which(x==-1)] = 2
  x[which(x==5)] = 3
  x[which(x==6)] = 4
  return(x)
}



data("sampleData")
# take a subset of the data for speed
seg <- CNSeg(sampleData[which(is.element(sampleData[, "ID"], 
                                         sample(unique(sampleData[, "ID"]), 10))), ])
rsBypair <- getRS(seg, by = "pair", imput = FALSE, XY = FALSE, what = "mean")
rsBypair
a = rs(rsBypair)

table(cnv$State)
class(a)
a[[1]]
dim(a[[1]])
head(a[[1]])
