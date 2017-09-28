library(qqman)

i=1

herdabilidades = c()

for(i in 1:22){
  
  cnv = read.table(paste("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/cromo", i, sep=""), header = T)
  cnv[1:10,1:10]
  
  
  a = read.table(paste("/home/cicconella/solar831/Chr", i, "/CNVs/resultados/herdabilidades",sep=""), sep="?")
  head(a)
  b = unlist(strsplit(as.character(a[,1]), split=":"))
  b = unlist(strsplit(b, split="p = "))
  b = unlist(strsplit(b, split="is "))
  b
  
  cns = b[seq(1,length(b),by=4)]
  herd = b[seq(3,length(b),by=4)]
  herd = as.numeric(herd)
  
  b = unlist(strsplit(b, split="(", fixed = T))
  
  
  pvalue = b[seq(4,length(b),by=5)]
  pvalue = as.numeric(pvalue)
  
  
  cns = unlist(strsplit(cns, split="CN"))
  cns = cns[seq(2,length(cns),by=2)]  
  
  final = as.data.frame(cbind(cns, pvalue, herd))
  head(final)
  final$cns = as.numeric(as.character(final$cns))
  head(final)
  final[final$cns==2,]
  final = final[order(final$cns),]
  head(final)
  
  plot(as.numeric(as.character(final$herd)))  
  
  head(final)
  
  final = cbind(cnv[,1:3], final)

  head(final)
    
  herdabilidades = rbind(herdabilidades, final)        
}

dim(herdabilidades)

head(herdabilidades)
colnames(herdabilidades) = c("Chr", "Start", "End","CNV", "P-value", "Herdabilidade")

class(herdabilidades)

head(herdabilidades)

herdabilidades$Herdabilidade = as.numeric(as.character(herdabilidades$Herdabilidade))

plot(herdabilidades$Herdabilidade, pch = 20, col="red")

h = herdabilidades[,c(1,2,6)]
colnames(h) = c("CHR", "BP", "P")
head(h)
h$P = as.numeric(as.character(h$P))
dim(h)

manhattan(h, col = c("darkgreen", "black"), main="CNVs Transmission Rate", 
          logp = F, ylim=c(0,1), ylab="Transmission Rate")

hist(h$P, col="darkgreen", nc=100)

a = length(which(h$P==0))
b = length(which(h$P==1))

p = h$P
p = p[-which(h$P==0)]
p
p = p[-which(p==1)]
p

his = hist(p, nc=20)
his$density = his$counts/sum(his$counts)*100
plot(his, freq=F, col="darkgreen", xlab = "Tansmission Rate", ylab = "Freq. (%)",
     main="Distribution of CNV Transmission Rate")

his$density

##### Novas herdabilidades #####

herdabilidades2 = c()
i=1
for(i in 1:22){
  
  cnv = read.table(paste("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/cromo", i, sep=""), header = T)
  cnv[1:10,1:10]
  
  a = read.table(paste("/home/cicconella/solar831/Chr", i, "/CNVs/resultados2/herdabilidades",sep=""), sep="?")
  head(a)
  dim(a)
  b = unlist(strsplit(as.character(a[,1]), split=":"))
  b = unlist(strsplit(b, split="p = "))
  b = unlist(strsplit(b, split="is "))
  b
  
  cns = b[seq(1,length(b),by=4)]
  herd = b[seq(3,length(b),by=4)]
  herd = as.numeric(herd)
  
  b = unlist(strsplit(b, split="(", fixed = T))
  
  
  pvalue = b[seq(4,length(b),by=5)]
  pvalue = as.numeric(pvalue)
  
  
  cns = unlist(strsplit(cns, split="CN"))
  cns = cns[seq(2,length(cns),by=2)]  
  
  final = as.data.frame(cbind(cns, pvalue, herd))
  head(final)
  final$cns = as.numeric(as.character(final$cns))
  head(final)
  final[final$cns==2,]
  final = final[order(final$cns),]
  head(final)
  
  final[final$cns==126,]
  
  plot(as.numeric(as.character(final$herd)))  
  
  head(final)
  
  final = cbind(cnv[,1:3], final)
  
  head(final)
  
  herdabilidades2 = rbind(herdabilidades2, final)        
}

dim(herdabilidades2)

head(herdabilidades2)
colnames(herdabilidades2) = c("Chr", "Start", "End","CNV", "P-value", "Herdabilidade")

class(herdabilidades2)

head(herdabilidades2)

herdabilidades2$Herdabilidade = as.numeric(as.character(herdabilidades2$Herdabilidade))

plot(herdabilidades2$Herdabilidade, pch = 20, col="red")

h = herdabilidades2[,c(1,2,6)]
colnames(h) = c("CHR", "BP", "P")
head(h)
h$P = as.numeric(as.character(h$P))
dim(h)

manhattan(h, col = c("darkgreen", "black"), main="CNVs Transmission Rate", 
          logp = F, ylim=c(0,1), ylab="Transmission Rate")

hist(h$P, col="darkgreen", nc=100)

a = length(which(h$P==0))
b = length(which(h$P==1))

p = h$P
p = p[-which(h$P==0)]
p
p = p[-which(p==1)]
p

his = hist(p, nc=20)
his$density = his$counts/sum(his$counts)*100
plot(his, freq=F, col="darkgreen", xlab = "Tansmission Rate", ylab = "Freq. (%)",
     main="Distribution of CNV Transmission Rate")

his$density

head(herdabilidades)
head(herdabilidades2)

table(herdabilidades$Herdabilidade==herdabilidades2$Herdabilidade)


a = herdabilidades2

head(a)

x = which(a$Start==72566893)

a[(x-9):(x+9),]
  

