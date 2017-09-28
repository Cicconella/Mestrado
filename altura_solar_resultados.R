library(qqman)

i=1

altura = c()

for(i in 1:22){
  
  cnv = read.table(paste("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/cromo", i, sep=""), header = T)
  cnv[1:10,1:10]
  
  
  a = read.table(paste("/home/cicconella/solar831/Chr", i, "/altura/resultados/p-values",sep=""), sep="?")
  head(a)
  b = unlist(strsplit(as.character(a[,1]), split=":"))
  b = unlist(strsplit(b, split="p = "))
  b = unlist(strsplit(b, split="(", fixed=T))
  b = b[seq(1,length(b),by=2)]
  cns = b[seq(1,length(b),by=2)]
  pvalue = b[seq(2,length(b),by=2)]

  cns = unlist(strsplit(cns, split="CN"))
  cns = cns[seq(2,length(cns),by=2)]
  final = as.data.frame(cbind(cns, pvalue))
  head(final)
  final$cns = as.numeric(final$cns)
  final = final[order(cns),]
  
  head(final)

  plot(final)  
  
  
  a = read.table(paste("/home/cicconella/solar831/Chr",i,"/altura/resultados/herdabilidades",sep=""), sep="?")
  head(a)
  b = unlist(strsplit(as.character(a[,1]), split=":"))
  b = unlist(strsplit(b, split="p = "))
  b = unlist(strsplit(b, split="is "))
  cns = b[seq(1,length(b),by=4)]
  herd = b[seq(3,length(b),by=4)]
  herd = as.numeric(herd)
  plot(herd)
  
  head(final) 

  final = cbind(cnv[,1:3], final, herd)

  head(final)
  
  altura = rbind(altura, final)        
}

dim(altura)

head(altura)
colnames(altura) = c("Chr", "Start", "End","CNV", "P-value", "Herdabilidade")

class(altura)

head(altura)
plot(altura$Herdabilidade, pch = 20, col="red")

altura[which(altura$Herdabilidade<0.8244920),]

altura$`P-value`= as.numeric(as.character(altura$`P-value`))

manAltura = altura[,c(4,1,2,5)]
colnames(manAltura) = c("SNP","CHR", "BP", "P")
head(manAltura)
manAltura$P = as.numeric(as.character(manAltura$P))
dim(manAltura)

manAltura$SNP = c(1:nrow(manAltura))

manhattan(manAltura, col = c("darkgreen", "black"), main="Height")
summary(manAltura)

altura[which((altura$`P-value`)==min(altura$`P-value`)),]
altura[which((altura$Herdabilidade)==min(altura$Herdabilidade)),]

head(altura)

manAltura = altura[,c(4,1,2,6)]
colnames(manAltura) = c("SNP","CHR", "BP", "P")
head(manAltura)
manAltura$P = as.numeric(as.character(manAltura$P))
dim(manAltura)

manhattan(manAltura, col = c("darkgreen", "black"), main="Height", logp = F, ylim=c(0.82,0.87), ylab="Heritability")
?manhattan
