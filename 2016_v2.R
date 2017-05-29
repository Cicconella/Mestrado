
########## Functions ##########

options(useHTTPS=FALSE, BioC_mirror="http://bioconductor.org")
source("http://bioconductor.org/biocLite.R")
biocLite("CNTools")
biocLite("Sushi")
library("CNTools")
library(Sushi)


itStwo = function(x){
  x[which(x==-1)] = 2 
  return(x)
}

cont_0 = function(x){
  length(which(x[-c(1,2,3)]==0))
}

cont_1 = function(x){
  length(which(x[-c(1,2,3)]==1))
}

cont_2 = function(x){
  length(which(x[-c(1,2,3)]==2))
}

cont_3 = function(x){
  length(which(x[-c(1,2,3)]==3))
}
cont_4 = function(x){
  length(which(x[-c(1,2,3)]==4))
}

mut_0 = apply(cnvs,1,cont_0)
mut_1 = apply(cnvs,1,cont_1)
mut_2 = apply(cnvs,1,cont_2)
mut_3 = apply(cnvs,1,cont_3)
mut_4 = apply(cnvs,1,cont_4)

mutacoes = cbind(mut_0,mut_1,mut_2,mut_3,mut_4)

grupos = function(a){
  t = sum(a)
  g = 5
  
  if(a[1]<t*maf)
    g = g-1
  if(a[2]<t*maf)
    g = g-1
  if(a[3]<t*maf)
    g = g-1
  if(a[4]<t*maf)
    g = g-1
  if(a[5]<t*maf)
    g = g-1
  
  return(g)
}


cont_0 = function(x){
  length(which(x[-c(1,2,3)]==0))
}
cont_1 = function(x){
  length(which(x[-c(1,2,3)]==1))
}

group = function(a){
  t = sum(a)
  g = 2
  
  if(a[1]<t*maf)
    g = g-1
  if(a[2]<t*maf)
    g = g-1
  
  return(g)
}


########## APT/PennCNV outputs ##########

### Information about CNV regions
cnv = read.table("../DadosMestrado/tableCNV", sep="\t", header = F)
colnames(cnv) = c("Chr","Start","End","Number","Length", "State", "CN", "Sample", "First Marker", "Last Marker")
head(cnv)
summary(cnv)
dim(cnv)
attach(cnv)

### Information about quality control PennCNV per sample
qc = read.table("../DadosMestrado/tableQC", sep="\t", header = T)
qc = qc[order(qc$File),] 
head(qc)
summary(qc)
dim(qc)

########## Cleaning Bad Samples ##########
length(which(qc$LRR_SD > 0.35))

length(which(qc$BAF_mean > 0.6))+length(which(qc$BAF_mean < 0.4))

length(which(qc$BAF_drift > 0.01))

length(which(qc$WF > 0.04))+length(which(qc$WF < -0.04))


fail1 = intersect(which(qc$LRR_SD > 0.35), which(qc$BAF_drift > 0.01))
fail2 = intersect(which(qc$LRR_SD > 0.35), which(qc$WF > 0.04))
fail3 = intersect(which(qc$LRR_SD > 0.35), which(qc$WF < 0.04))

head(cnv)
dim(cnv)

qc = qc[-c(fail2,fail3),]
dim(qc)
head(qc)
attach(qc)

cnv = merge(cnv, qc, by.x = "Sample", by.y = "File")
cnv = cnv[,-c(11:19)]
dim(cnv)
head(cnv)
summary(cnv)
attach(cnv)
ind = unique(Sample)
length(ind)

remove(qc)
remove(fail3)
remove(fail2)
remove(fail1)
#remove(ind)

########## BY chromosome ##########

# Find the minimal regions for each chromosome

chromosomes <- vector("list",22)

for(i in 1:22){
  
  cnvA = cnv[Chr==i,]
  dim(cnvA)
  head(cnvA)
  aux = cnvA[order(cnvA$Start),]
  dim(aux)
  head(aux)
  
  # Change matrix to data.frame and add an auxiliar region from the start of the 
  # chromossome
  
  aux <- data.frame(c(aux$Sample,ind), c(aux$Chr,rep(i, length(ind))), 
                    c(aux$Start,rep(1, length(ind))), c(aux$End,rep(max(aux$End), length(ind))), 
                    c(aux$Number,rep(1, length(ind))), c(aux$CN,rep(-1, length(ind))))
  colnames(aux) = c("ID", "chrom", "loc.start", "loc.end", "num.mark",  "seg.mean")
  str(aux)
  head(aux)
  dim(aux)
  
  seg <-  CNSeg(aux)
  seg
  rsByregion <- getRS(seg, by = "region", imput = TRUE, XY = FALSE, what = "max")
  cnvA = rs(rsByregion)
  head(cnvA)
  
  dim(cnvA)
  cnvA = apply(cnvA,2,itStwo)
  
  chromosomes[[i]] = cnvA
}

remove(cnv)
remove(cnvA)
remove(aux)
remove(rsByregion)
remove(seg)

chromosomes[[1]][1,]


contaMut = function(x){
  return(length(which(x!=2))-3)
}

contaMut(chromosomes[[1]][2,])

plotar = cbind(chromosomes[[1]][,1:3], apply(chromosomes[[1]], 1, contaMut))
plotar = as.data.frame(plotar)
colnames(plotar) = c("chrom", "start", "end", "score")
plotar = plotar[-1,]
plotar[,1] = "chr1"

plotar[1:10,]
tail(plotar)

apply(plotar, 2, class)

chrom = plotar[1,1]
chromstart = chromosomes[[1]][2,2]
chromend = chromosomes[[1]][nrow(chromosomes[[1]]),3]



plotar[which(plotar[,4]>800),]
plotar[which(plotar[,2]>110000000),]

png("../DadosMestrado/Sushi/chr1.png")


plotBedgraph(plotar, chrom, chromstart,
             chromend, main = "CNVs per region (1019 samples)",
             colorbycol= SushiColors(5))

labelgenome(chrom, chromstart,chromend,n=4,scale="Mb")
mtext("Number of CNVs",side=2,line=2.5,cex=1,font=2)
axis(side=2,las=2,tcl=.2)

dev.off()

for(i in 1:22){
  plotar = cbind(chromosomes[[i]][,1:3], apply(chromosomes[[i]], 1, contaMut))
  
  plotar = as.data.frame(plotar)
  colnames(plotar) = c("chrom", "start", "end", "score")
  
  plotar[,1] = paste("chr",i, sep="")
  
  plotar$score[which(plotar$score<0)]=0
  
  chrom = plotar[1,1]
  chromstart = chromosomes[[i]][2,2]
  chromend = chromosomes[[i]][nrow(chromosomes[[i]]),3]
  
  png(paste("../DadosMestrado/Sushi/chr",i,".png", sep=""))
      plotBedgraph(plotar, chrom, chromstart,
                   chromend, main = "CNVs per region (1019 samples)",
                   colorbycol= SushiColors(5))
      
      labelgenome(chrom, chromstart,chromend,n=4,scale="Mb")
      mtext("Number of CNVs",side=2,line=2.5,cex=1,font=2)
      axis(side=2,las=2,tcl=.2)
  dev.off()
  
  print(i)
}





sizesT= c()
for(x in 1:22){
  y = dim(chromosomes[[x]])[1]
  sizesT = c(sizesT,y)
}

png("/Users/Ana/Google Drive/2016/before.png")
plot(1:22,sizesT, xlab="Chromosomes", ylab="Number of CNV Regions", 
     main = "CNV Regions before Cleaning",pch=16)
dev.off()

remove(x)
remove(y)
remove(i)

########## Cleaning CNV Regions ##########

for(i in 1:22){
  mut_0 = apply(chromosomes[[i]],1,cont_0)
  mut_1 = apply(chromosomes[[i]],1,cont_1)
  mut_2 = apply(chromosomes[[i]],1,cont_2)
  mut_3 = apply(chromosomes[[i]],1,cont_3)
  mut_4 = apply(chromosomes[[i]],1,cont_4)
  
  mutations = cbind(mut_0,mut_1,mut_2,mut_3,mut_4)
  
  dim(chromosomes[[i]])
  dim(mutations)
  
  groupCNV = apply(mutations, 1, grupos)
  
  summary(as.factor(groupCNV))
  
  exclude = which(groupCNV==1)
  chromosomes[[i]] = chromosomes[[i]][-exclude,]
  dim(chromosomes[[i]])
  head(chromosomes[[i]][,1:15])
  print(i)
}

sizes = c()
for(x in 1:22){
  y = dim(chromosomes[[x]])[1]
  sizes = c(sizes,y)
}


########## Cleaning CNV Regions - Binary##########

maf = 0.05

chrBin = vector("list",22)

for(i in 1:22){
  chrBin[[i]] = chromosomes[[i]]
  chrBin[[i]][which(chrBin[[i]]!=2)] = 1
  chrBin[[i]][which(chrBin[[i]]==2)] = 0
  
  chrBin[[i]][,c(1,2,3)] = chromosomes[[i]][,c(1,2,3)]
}

head(chromosomes[[1]][,1:10])
head(chrBin[[1]][,1:10])


for(i in 1:22){
  mut_0 = apply(chrBin[[i]],1,cont_0)
  mut_1 = apply(chrBin[[i]],1,cont_1)
  
  mutations = cbind(mut_0,mut_1)
  
  dim(chrBin[[i]])
  dim(mutations)
  
  groupCNV = apply(mutations, 1, group)
  
  summary(as.factor(groupCNV))
  
  exclude = which(groupCNV==1)
  chrBin[[i]] = chrBin[[i]][-exclude,]
  dim(chrBin[[i]])
  head(chrBin[[i]][,1:15])
  print(i)
}

sizesB = c()
for(x in 1:22){
  y = dim(chrBin[[x]])[1]
  sizesB = c(sizesB,y)
}

png("/Users/Ana/Google Drive/2016/after-bin-maf-0.02.png")
plot(1:22,sizesB, xlab="Chromosomes", ylab="Number of CNV Regions", 
     main = "CNV Regions after Cleaning",pch=16)
dev.off()

png("/Users/Ana/Google Drive/2016/befofe-after-bin-0.02.png")
plot(sizesB,sizesT, xlab="Number of Filtered Regions", ylab = "Total of Regions",
     xlim=c(-20,2400),ylim=c(-20,11000),pch=16)
identify(sizesB,sizesT)
dev.off()

remove(mutations)
remove(exclude)
remove(groupCNV)
remove(i)
remove(maf)
remove(mut_0)
remove(mut_1)
remove(x)
remove(y)plotar[,2]>125000000

########## Individuals ##########

# Sample Data

info = read.table("/home/cicconella/Dropbox/2016/Project/dados2", header = T, sep = ",")

head(info)  
summary(info)
dim(info)
str(info)
attach(info)

# Association between sample id and celfiles

ind = read.table("/home/cicconella/Dropbox/2016/Project/individuos")
colnames(ind) = c("cel", "IID")

class(ind)
head(ind)
dim(ind)
summary(ind)
str(ind)

# Associacao 

cel = IID Vicious Delicious

for(i in 1:length(cel)){
  
  aux = which(ind$IID==cel[i])
  if(length(aux)!=0)
    cel[i] = max(ind$cel[aux])
  else
    cel[i] = NAn = length(unique(dataset$FID))
}

cel[1:10]

#  Info + cel

info = cbind(info,cel)

head(info)  
summary(info)
dim(info)
str(info)


# Correct CEL Files
info[which(info$IID==4919),9]=1444
info[which(info$IID==15908),9]=1455
info[which(info$IID==15911),9]=2113
info[which(info$IID==15921),9]=482
info[which(info$IID==32608),9]=1457
info[which(info$IID==96502),9]=2360

attach(info)

head(info)

positions = rep(0,nrow(info))
names = as.numeric(colnames(chromosomes[[1]])[-c(1,2,3)])

for (i in 1:length(IID)){
  k = which(info$cel[i]==names)
  if(length(k)==1)
    positions[i] = k
  else
    positions[i] = NA
}

remove(ind)
remove(aux)
remove(i)
remove(names)
remove(k)
remove(cel)

########## Getting the files ped and phen ##########

# .ped

ped = cbind(info$IID,info$PAT,info$MAT,info$SEX,info$FID)
colnames(ped) = c("id","fa","mo","sex","fid")

write.table(ped, "/Users/Ana/Google Drive/2016/files/samples.ped", row.names = F, quote = F, sep = ",")

# .phen

phen = vector("list",22)

for(j in 1:22){
  fen = c()
  for(i in 1:nrow(chrBin[[j]])){  
    print(c("CHR ", j, "Reg ",i))
    cn = chrBin[[j]][i,-c(1,2,3)]
    cn = cn[positions]
    
    fen = cbind(fen,cn)
  }
  phen[[j]] = cbind(info,fen)
}

dim(phen[[1]])
dim(chrBin[[1]])

for(i in 1:22){
  
  phenotypes = phen[[i]][,-c(1,3,4,6,9)]
  
  name = paste("/Users/Ana/Google Drive/2016/files/phen",i,".phen",sep="")
  write.table(phenotypes,name,row.names = F, quote = F, sep = ",")
}

remove(fen)
remove(info)
remove(ped)
remove(phenotypes)
remove(chrBin) #ja esta organizado nos arquivos .phen
remove(cn)
remove(i)
remove(j)
remove(name)
remove(positions)
remove(sizesB)
remove(sizesT)
remove(phen)

########## Correctins files to Solar##########

for(i in 2:22){
  
  name = paste("/Users/Ana/Google Drive/2016/files/phen",i,".phen",sep="")
  chr = read.table(name, header = T, sep=",")
  
  chr[1:10,1:10]
  dim(chr)
  
  colnames(chr)[1:3] = c("id", "sexo", "idade")
  
  for(i in 1:nrow(chr)){
    chr[i,which(is.na(chr[i,]))] = ""
  }
  
  write.table(chr, name, row.names = F, quote = F, sep = ",")
}


rm(name)

