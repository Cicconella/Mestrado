
########## Functions ##########

options(useHTTPS=FALSE, BioC_mirror="http://bioconductor.org")
source("http://bioconductor.org/biocLite.R")
#biocLite("CNTools")
#biocLite("Sushi")
library("CNTools")
library(Sushi)
library(stargazer)
library("kinship2")

arruma = function(x){
  x[which(x==1)] = 0
  x[which(x==2)] = 1
  x[which(x==-1)] = 2
  x[which(x==5)] = 3
  x[which(x==6)] = 4
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



group = function(a){
  t = sum(a)
  g = 2
  
  if(a[1]<t*maf)
    g = g-1
  if(a[2]<t*maf)
    g = g-1
  
  return(g)
}


contaMut = function(x){
  return(length(which(x!=2))-3)
}


getwd()


##### Dados CNV #####

dir = "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado"

a = read.table(paste(dir,"/tableCNV", sep=""))
head(a)

colnames(a) = c("Chr", "Begin", "End", "NumberSNPs", "Size", "State", "CN",
                "Sample", "FirstMarker", "LastMarker")
head(a)
dim(a)
length(unique(a$Sample))
sum(table(a$Sample))

##### Individuos Limpos #####

ind = read.table(paste(dir,"/ind_limpos", sep=""), header=T)
head(ind)

individuos = ind$cel
individuos = as.data.frame(individuos)
head(individuos)
head(a)
dim(a)

a = merge(a, individuos, by.x = "Sample", by.y = "individuos")
head(a)
dim(a)

attach(a)
head(a)


##### Samples #####
Sample
print(paste("Total de amostras:", length(unique(Sample))))

summary(as.numeric(table(Sample)))
sd(as.numeric(table(Sample)))

#png(paste(dir, "/numberCNVs.png", sep=""))
hist(table(Sample), nc = 1000, main = "Number of CNVs per sample", xlab = "Number of CNVs")
#dev.off()

table(Sample)[which(table(Sample)>500)]
paste(length(which(table(Sample)>500))/length(unique(Sample)), "são maiores que 500")
paste(length(which(table(Sample)>250))/length(unique(Sample)), "são maiores que 250")
paste(length(which(table(Sample)>100))/length(unique(Sample)), "são maiores que 100")

table(Sample)
summary(as.numeric(table(Sample)))

normal = table(Sample)
length(normal)
normal = normal[normal<101]
length(normal)
753/910
summary(as.numeric(normal))
sd(as.numeric(normal))
hist(as.numeric(normal), nc=100)
shapiro.test(as.numeric(normal))
shapiro.test(rnorm(100,5,2))

porcentagens = c()

for(i in seq(0,3000, by =25)){
  porcentagens = c(porcentagens, length(which(table(Sample)<i))/length(unique(Sample)))
}

porcentagens

#plot(seq(0,4000, by =25), porcentagens, pch="", ylab="% Samples", xlab = "Number of CNVs")
#png(paste(dir, "/samplesSize.png", sep=""))
#lines(seq(0,4000, by =25), porcentagens, main = "Frequency of samples by number of CNVs")
#dev.off()

bla = cbind(seq(0,3000, by =25), porcentagens)

print(paste("% de pessoas com menos de 100 CNVs:", porcentagens[5]))

stargazer(bla, digits = 2)

head(a)

#png(paste(dir, "/absCNV.png", sep=""))
barplot(table(CN)/sum(table(CN)), col = "blue", xlab = "Copy Number", 
        ylab = "Frequency (%)", main = "Frequency of Copy Number",
        ylim = c(0,0.7))
#dev.off()

table(CN)/sum(table(CN))

which(table(Sample)<101)
names(which(table(Sample)<101))

boxplot(as.numeric(table(Sample)[table(Sample)<101]), col = "blue",
        main="Distribution of Number of CNVs per Sample", pch=20)

bla = a[a$Sample %in% names(which(table(Sample)<101)),]

bla[1:100,]

png(paste(dir, "/absCNV100.png", sep=""))
barplot(table(bla$CN)/sum(table(bla$CN)), col = "blue", xlab = "Copy Number", 
        ylab = "Frequency (%)", main = "Frequency of Copy Number",
        ylim = c(0,0.7))
dev.off()

table(bla$CN)/sum(table(bla$CN))


head(a)
length(unique(a$Sample))

infos = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/informacoesIndividuos", header=T)
head(infos)

head(a)

cnv_idade = merge(a,infos, by.x="Sample", by.y="cel")

head(cnv_idade)
length(unique(cnv_idade$Sample))

cnv_idade=cnv_idade[,c("Sample","Idade")]
head(cnv_idade)
length(unique(cnv_idade$Sample))

head(cnv_idade)
tab=table(cnv_idade$Sample)
base = unique(cnv_idade)
head(base)
cnv_idade = cbind(names(tab),as.numeric(tab),base)
head(cnv_idade)
plot(cnv_idade$Idade,cnv_idade$`as.numeric(tab)`, pch=20)

cnv_idade = cnv_idade[cnv_idade$`as.numeric(tab)`<100,]
head(cnv_idade)
linear=lm(cnv_idade$`as.numeric(tab)`~cnv_idade$Idade)
linear
plot(cnv_idade$Idade,cnv_idade$`as.numeric(tab)`, pch = 20, col="blue",
     xlab="Age", ylab="Number of CNVs", main="Number of CNVs per Age")
abline(linear)
anova(linear)

##### Size #####
head(a)

summary(Size)
log10(3)
log10(27435314)

hist(Size, nc=100)
barplot(table(Size))

mean(Size)
summary(Size)

data = Size

log10(summary(Size))
summary(log10(data))
10**summary(log10(data))
10**(1:10)

#png(paste(dir, "/totalLength.png", sep=""))
hist(log10(data), nc=100, xlab = "Length of CNV", ylab = "Absolute Frequency", col="blue", main = "Histogram of CNV Length", ylim = c(0,6000))
hist(log10(data), nc=100,xaxt="n", xlab = "Length of CNV", ylab = "Absolute Frequency", col="blue", main = "Histogram of CNV Length", ylim = c(0,6000))
axis(1, at=1:7, labels = c("10bp", "100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb"))
#dev.off()

head(a)

#png(paste(dir, "/totalDel.png", sep=""))
hist(log10(Size[a$CN<2]), 
     nc=100,xaxt="n", xlab = "Length of CNV", ylab = "Absolute Frequency", 
     col="blue", main = "Histogram of Deletions Length", ylim = c(0,6000))
axis(1, at=1:7, labels = c("10bp", "100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb"))
#dev.off()

#png(paste(dir, "/totalDup.png", sep=""))
hist(log10(Size[a$CN>2]), 
     nc=100,xaxt="n", xlab = "Length of CNV", ylab = "Absolute Frequency", 
     col="blue", main = "Histogram of Duplications Length", ylim = c(0,6000))
axis(1, at=1:7, labels = c("10bp", "100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb"))
#dev.off()


########## BY chromosome ##########

head(a)
head(ind)

# Find the minimal regions for each chromosome

chromosomes <- vector("list",22)

i=22

for(i in 1:22){
  
  cnvA = a[Chr==i,]
  dim(cnvA)
  head(cnvA)
  aux = cnvA[order(cnvA$Begin),]
  dim(aux)
  head(aux)
  
  # Change matrix to data.frame and add an auxiliar region from the start of the 
  # chromossome
  
  aux <- data.frame(c(aux$Sample,unique(Sample)), c(aux$Chr,rep(i, nrow(ind))), 
                    c(aux$Begin,rep(1, nrow(ind))), c(aux$End,rep(max(aux$End), nrow(ind))), 
                    c(aux$Number,rep(1, nrow(ind))), c(aux$State,rep(-1, nrow(ind))))
  colnames(aux) = c("ID", "chrom", "loc.start", "loc.end", "num.mark",  "seg.mean")
  str(aux)
  head(aux)
  dim(aux)
  
  aux[aux$ID==2,]
  
  
  seg <-  CNSeg(aux)
  seg
  rsByregion <- getRS(seg, by = "region", imput = TRUE, XY = FALSE, what = "max")
  cnvA = rs(rsByregion)
  cnvA[1:10,1:10]
  
  dim(cnvA)
  cnvA = cbind(cnvA[,1:3],apply(cnvA[,-c(1:3)],2,arruma))
  
  cnvA[200:400,1:5]
  
  
  
  chromosomes[[i]] = cnvA
  
  
}

remove(cnvA)
remove(aux)
remove(rsByregion)
remove(seg)

total = c()

for(i in 1:22){
 total = c(total, dim(chromosomes[[i]])[1])  
}

total
sum(total)



barplot(table(as.numeric(chromosomes[[1]][100,-c(1:3)])), col="blue", ylab="Number of Samples",
        xlab="CN Group", main = "Number of Samples x Number of Copies")

910*0.02

chromosomes[[1]][100,1:10]

head(a)

dim(a)



##### CNVs, onde estao? #####

table(a$Chr)
barplot(table(a$Chr), xlab = "Chromosome", names.arg=c(1:22),col="blue", cex.names  = 0.75,
        main = "CNVs detected by Chromosome", ylab = "Absolute frequency")


tamanhos = read.table(paste(dir, "/tamanhos", sep=""), header = F, sep = "\t")
head(tamanhos)

prop = c()
prop2 = c()

i=1

for(i in 1:22){
  
  if(i!=2){
    plotar = cbind(chromosomes[[i]][,1:3], apply(chromosomes[[i]], 1, contaMut))  
  } else{
    plotar = cbind(chromosomes[[i]][,1:3], (apply(chromosomes[[i]], 1, contaMut)+1))    
  }
  
  
  plotar = as.data.frame(plotar)
  colnames(plotar) = c("chrom", "start", "end", "score")
  
  plotar[,1] = paste("chr",i, sep="")
  
  chrom = plotar[1,1]
  chromstart = chromosomes[[i]][1,2]
  chromend = chromosomes[[i]][nrow(chromosomes[[i]]),3]
  
  # png(paste(dir, "/Sushi/chr",i,".png", sep=""))
  # plotBedgraph(plotar, chrom, chromstart,
  #              chromend, main = "CNVs per region (910 samples)",
  #              colorbycol= SushiColors(5), range = c(0,910))
  # 
  # labelgenome(chrom, chromstart,chromend,n=4,scale="Mb")
  # mtext("Number of CNVs",side=2,line=2.5,cex=1,font=2)
  # axis(side=2,las=2,tcl=.2)
  # dev.off()
  
  print(i)
  
  
  bla = plotar$end-plotar$start+1
  bla = sum(bla[which(plotar$score!=0)])
  
  prop = c(prop, bla/tamanhos$V2[i])
  
  bla = plotar$end-plotar$start+1
  bla = sum(bla[which(plotar$score>20)])
  
  prop2 = c(prop2, bla/tamanhos$V2[i])
  
}

plotar
plotar[which(plotar$score==max(plotar$score)),]

##### Teste #####

# head(plotar)
# 
# x=which(plotar$score==max(plotar$score))
# 
# plotar[(x-20):(x+20),]
# plotar[(x-12):(x+10),]
# mean(plotar[(x-12):(x+10),4])
# 
# chromosomes[[i]][x,1:3]
# contaMut(chromosomes[[i]][x,1:3])
# 
# stargazer(as.numeric(table(as.numeric(chromosomes[[i]][x,-(1:3)]))))
# 
# teste = 764792
# teste2 =802515
# plotar$start-764792
# 
# plotar[20,]


sizesT= c()
for(x in 1:22){
  y = dim(chromosomes[[x]])[1]
  sizesT = c(sizesT,y)
}

png(paste(dir,"/plots/before.png", sep=""))
plot(1:22,sizesT, xlab="Chromosomes", ylab="Number of CNV Regions", 
     main = "CNV Regions before Cleaning",pch=16)
dev.off()

########## Cleaning CNV Regions ##########

i=1

original = chromosomes

dim(original[[1]])

#chromosomes = original
original[[1]][1:10,1:10]
chromosomes[[1]][1:10,1:10]

maf = 0.02

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
for(i in 1:22){
  y = dim(chromosomes[[i]])[1]
  sizes = c(sizes,y)
}

sum(sizes)

barplot(table(as.numeric(chromosomes[[1]][100,-c(1:3)])), col="blue", ylab="Number of Samples",
        xlab="CN Group", main = "Number of Samples x Number of Copies")

910*0.02

chromosomes[[1]][100,1:10]


plot(sizes/total)
plot(total)


png(paste(dir, "/plots/before2.png", sep=""))
plot(1:22,sizes, xlab="Chromosomes", ylab="Number of CNV Regions", 
     main = "CNV Regions after Cleaning",pch=16)
dev.off()


bla = rbind(cbind(1:22,sizes),cbind(1:22,sizesT))

png(paste(dir, "/plots/CNVnumber.png", sep=""))
plot(bla, xlab="Chromosomes", ylab="Number of CNVs", 
     main = "Number of CNVs", pch=16, col = c(rep("blue", 22),rep("red", 22)))
legend("topright", c("Before cleaning", "After cleaning"), col = c("red", "blue"), pch=16)
dev.off()


for(i in 1:22)
  write.table(chromosomes[[i]], paste(dir, "/Cromossomos/cromo",i, sep=""), row.names = F)



##### 1000Genomes


c1 = read.table(paste(dir, "/1000genomes/c1", sep=""), header = T)

head(c1)

cnv = read.table(paste("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/cromo", 1, sep=""), header = T)
cnv[1:10,1:10]

head(c1)







##### Tamanho de CNVs Minimais #####

total  = c()

head(chromosomes[[1]][1:10,1:10])

for(i in 1:22){
  total = c(total, chromosomes[[i]][,3]-chromosomes[[i]][,2])  
}
total= total+1

length(total)

summary(total)

hist(log10(total))
hist(log10(total), nc=100, xlab = "Length of CNV", 
     ylab = "Absolute Frequency", col="blue", main = "Histogram of CNV Length", 
     ylim = c(0,1000))
hist(log10(total), nc=100,xaxt="n", xlab = "Length of CNV", 
     ylab = "Absolute Frequency", col="blue", main = "Histogram of CNV Length", 
     ylim = c(0,1000))
axis(1, at=1:7, labels = c("10bp", "100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb"))




##### Quantidade por individuo #####

chromosomes[[1]][1:10,1:10]


contaMut = function(x){
  return(length(which(x!=2)))  
}

numero = c()

for(i in 1:22){
  numero=rbind(numero,apply(chromosomes[[i]][,-c(1:3)],2,contaMut))
}

dim(numero)
numero[1:10,1:10]

numero = apply(numero,2,sum)
numero
hist(numero, nc=100, col="blue", main="Number of CNVs per Sample", xlab="Number of CNVs")
summary(numero)

length(numero[numero<1001])
801/910

##### Genotipo e fenotipo #####

chromosomes[[1]][1:10,1:10]
table(as.numeric(chromosomes[[1]][4,-c(1:3)]))

ind = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/informacoesIndividuos", header=T)

head(ind)

ind[which(ind$cel==1),]

chromosomes[[1]][1:10,1:10]

head(ind)

chromosomes[[1]][1,1:50]

dim(chromosomes[[4]])

seq(1,676, by=10)

## Seleciona o genotipo

somas = c()
medias = c()

k=22

for(k in 1:22){
  
  sequencia = nrow(chromosomes[[k]])
  
  sequencia = sample(seq(1,sequencia))
  
  
  for(c in sequencia){
    genotipo = as.numeric(chromosomes[[k]][c,-c(1:3)])
    
    genotipo = (cbind(colnames(chromosomes[[k]])[-c(1:3)], genotipo))
    
    colnames(genotipo) = c("celfiles", "cnv")
    genotipo = as.data.frame(genotipo)
    head(genotipo)
    
    head(ind)
    
    final = merge(ind, genotipo, by.x = "cel", by.y = "celfiles", all.x = T)
    final[1:100,]
    tail(final)
    
    head(final)
    dim(final)
    
    comb = expand.grid(c(0:4),c(0:4),c(0:4))
    comb = data.frame(cbind(comb, rep(0,125)))
    comb[,4] = as.numeric(as.character(comb[,4]))
    
    colnames(comb) = c("P1", "P2", "OF", "CN")
    
    head(comb)
    head(final)
    
    i=1
    
    for(i in 1:nrow(final)){
      final[i,]
      
      if(is.na(as.character(final[i,"cnv"]))){
        #print("next")
        next
      }
      
      if(final[i,4]!=0 & final[i,5]!=0){
        if(is.na(as.character(final[which(final[,3]==final[i,4]),"cnv"]))
           | is.na(as.character(final[which(final[,3]==final[i,5]),"cnv"]))){
          next
        }else{
          a = t(as.matrix(c(as.numeric(as.character(final[which(final[,3]==final[i,4]),10])),
                            as.numeric(as.character(final[which(final[,3]==final[i,5]),10])),
                            as.numeric(as.character(final[i,10])))))
          colnames(a) = colnames(comb)[-4]
          comb[which(apply(comb, 1, function(x) identical(x[1:3], a[1,]))),4] = comb[which(apply(comb, 1, function(x) identical(x[1:3], a[1,]))),4]+1
        }
      }else{
        #print("orfao")
      }
    }
    
    #png(paste(getwd(),"/trios.png", sep = ""), width = 1400, height = 460)
    #barplot(comb[,4],pch=16, names.arg = apply(comb[,1:3],1,paste,collapse = ""), las=2)
    #dev.off()
    
    for(i in 1:75){
      a = t(as.matrix(as.numeric(c(rev(comb[i, 1:2]),comb[i,3]))))
      colnames(a) = c("P1", "P2","OF")
      bla = which(apply(comb, 1, function(x) identical(x[1:3], a[1,])))
      if(bla!=i){
        comb[i,4] = comb[i,4]+comb[bla,4]
        comb = comb[-bla,]
      }
    }
    
    dim(comb)
    
    #png(paste(getwd(),"/trios2.png", sep = ""), width = 1400, height = 460)
    #barplot(comb[,4],pch=16, names.arg = apply(comb[,1:3],1,paste,collapse = ""), las=2)
    #dev.off()
    
    head(comb)
    
    if(c == sequencia[1]){
      trios = comb
    }else{
      trios = cbind(trios,comb[,4])  
    }
    print(c)
    
  }
  head(trios)
  trios[,1:20]
  dim(trios)
  
  trios[,-c(1:3)]=trios[,-c(1:3)]/unique(apply(trios[,-(1:3)],2,sum))
    
  combinacoes = apply(trios[,-(1:3)], 1,mean)
  names(combinacoes) = apply(trios[,1:3],1,paste,collapse="")
  medias = cbind(medias,combinacoes)
  combinacoes_limpo = combinacoes[-which(names(combinacoes)=="222")] 
  combinacoes*100

    
  png(paste0("mediaTriosChr",k,".png"), width = 1200, height = 540)
  barplot(combinacoes_limpo, las=2, xlab = "Genotype", ylab = "Mean occurance per CNV",
          col = "blue", main = paste0("Occurance of CNVs in Trios (Chromosome ",k,")"),
          ylim = c(0,1))
  dev.off()
    
  combinacoes = apply(trios[,-(1:3)], 1,sum)
  names(combinacoes) = apply(trios[,1:3],1,paste,collapse="")
  somas = cbind(somas,combinacoes)
  combinacoes_limpo = combinacoes[-which(names(combinacoes)=="222")] 
  combinacoes

  
      
  png(paste0("somaTriosChr",k,".png"), width = 1200, height = 540)
  barplot(combinacoes_limpo, las=2, xlab = "Genotype", ylab = "Total occurance",
          col = "blue", main = paste0("Occurance of CNVs in Trios (Chromosome ",k,")"),
          ylim = c(0,200))
  dev.off()
    
}

#aux = medias
dim(medias)
head(medias)
head(somas)
colnames(somas) = paste("Chr.", seq(1,22))

head(somas)


medias = cbind(medias, apply(medias,1,mean),apply(medias,1,sd))
#colnames(medias) = c(paste("Chr.", seq(1,22)),"Mean")

barplot(apply(medias,1,mean))

medias = medias*100
colnames(medias) = c(paste("Chr.", seq(1,22)), "Mean", "Std. dev.")

head(medias)
stargazer(t(medias)[,1:15], summary = F, digits = 2)
stargazer(t(medias)[,16:30], summary = F, digits = 2)
stargazer(t(medias)[,31:45], summary = F, digits = 2)
stargazer(t(medias)[,46:60], summary = F, digits = 2)
stargazer(t(medias)[,61:75], summary = F, digits = 2)

dim(trios)

trios[,1:20]

trios[,1:15]

x = cbind(apply(trios[,1:3],1,paste,collapse=""),trios[,4:15])
colnames(x) = c("CNs", colnames(x)[-1])

stargazer(x, summary = F)
