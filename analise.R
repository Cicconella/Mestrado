########## Selecionar resultados ##########

# Qual MAF? 2: 0.02 e 5: 0.05
maf = 5

final = c()

for(i in 1:22){
  nome = paste("/home/cicconella/Dropbox/2016/binaryAnalysis/Saidas",maf,"/saida",i,"/herPvalue.txt", sep="")
  
  aux = read.table(nome, sep = "%")
  aux = as.matrix(aux)
  
  aux = unlist(strsplit(aux, split = "H2r is"))
  aux = unlist(strsplit(aux, split = "p ="))
  aux = unlist(strsplit(aux, split = "(Not", fixed = TRUE))
  aux = unlist(strsplit(aux, split = "(S", fixed = TRUE))
  
  aux[1:10]
  
  aux = matrix(aux, nrow = 4)
  
  h = as.numeric(aux[2,])
  p = as.numeric(aux[3,])
  p = -log(p, base=10)
  
  nome = paste("/home/cicconella/Dropbox/2016/binaryAnalysis/Saidas",maf,"/saida",i,"/trait.txt", sep="")
  aux = read.table(nome, sep = "%")
  aux = as.matrix(aux)
  aux[which(aux=="cn")]="cn.0"
  
  aux = unlist(strsplit(aux, split = "cn."))
  trait = as.numeric(aux[seq(2,length(aux),by=2)])
  
  chr = rep(i,length(trait))
  
  resultados = cbind(chr, trait, h, p)
  resultados = resultados[order(resultados[,2]),]
  
  final = rbind(final, resultados)
  
}

final = as.data.frame(final)

summary(final) 
head(final)
attach(final)

plot(final$h, final$p, xlab = "Heritability", ylab = "-log10(P-value)", col="blue", pch =16, main = "Heritability x -log10(P-Value)")
hist(h)

nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados CN/man","pCN",sep="")

jpeg(nome)
manH = as.data.frame(cbind(CN,Chr,pos,pCN))
colnames(manH)=c("SNP", "CHR", "BP","P")  
head(manH)
manhattan(manH, col = c("blue", "orange"))
dev.off()


########## Selecionar resultados ##########

final0 = final[which(final$h==0),] 
final05 = final[which(final$h>0.3 & final$h<0.5),] 
final1 = final[which(final$h==1),] 

head(final0)
dim(final0)
head(final1)
head(final05)

# Selecionar 3 CNVs de cada grupo

samples0 = sort(round(runif(3, 0, nrow(final0))))
a = final0[samples0,]

samples05 = sort(round(runif(3, 0, nrow(final05))))
b = final05[samples05,]

samples1 = sort(round(runif(3, 0, nrow(final1))))
c = final1[samples1,]

d = list(a,b,c)

# Para cada CNV selecionado, obter os fenotipos, associar as familias e plotar a proporçao de pessoas afetadas.

dataset = read.table("/home/cicconella/Dropbox/2016/Project/dataset", header = T)

for(j in 1:3){
  nome = paste("/home/cicconella/Dropbox/2016/binaryAnalysis/ex",j,".jpg", sep="")
  jpeg(nome,width = 800, height = 300)
  par(mfrow=c(1,3))
    for(i in 1:3){
    nome = paste("/home/cicconella/Dropbox/2016/files/",maf,"phen",d[[j]][i,1],".phen", sep="")
    aux = read.table(nome, sep=",", header = T)
    cn = aux[,5+d[[j]][i,2]]
    cn[which(is.na(cn))] = 2
    summary(cn)
    cn = as.factor(cn)
    info = as.data.frame(cbind(dataset[,1],cn))
    
    # 1=Normal 2=Affected 3=NA
    tabela = table(info)
    head(tabela)
    
    total = tabela[,1]+tabela[,2] 
    head(total)
    
    prop = tabela[,2]/total
    prop[which(is.na(prop))] = 0
    nome = paste("CNV from chromosome", d[[j]][i,1])
    plot(prop, xlab = "Family ID", ylab = "Proportion", pch=16, col="blue", main=nome)
    lines(seq(0,122,by=0.1),rep(mean(prop),length(seq(0,122,by=0.1))), col="red")
  }
  dev.off()
}

rm(list = ls())

########## Selecionar resultados da regressao altura covariantes ##########
par(mfrow=c(1,1))

final = c()

for(i in 1:22){
  nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados Covariaveis/resultado", i,"/copy", sep="")
  
  aux = read.table(nome, sep = "%")
  aux = as.matrix(aux)
  
  aux = unlist(strsplit(aux, split = "cn."))
  aux = unlist(strsplit(aux, split = "p ="))
  aux = unlist(strsplit(aux, split = "(Not", fixed = TRUE))
  aux = unlist(strsplit(aux, split = "(S", fixed = TRUE))
  
  aux[1:10]
  
  aux = matrix(aux, nrow = 4)
  
  cn = as.numeric(aux[2,])
  pcn = as.numeric(aux[3,])
  
  nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados Covariaveis/resultado", i,"/hp", sep="")
  
  aux = read.table(nome, sep = "%")
  aux = as.matrix(aux)
  
  aux = unlist(strsplit(aux, split = "H2r is"))
  aux = unlist(strsplit(aux, split = "p ="))
  aux = unlist(strsplit(aux, split = "(Not", fixed = TRUE))
  aux = unlist(strsplit(aux, split = "(S", fixed = TRUE))
  
  aux[1:10]
  
  aux = matrix(aux, nrow = 4)
  
  h = as.numeric(aux[2,])
  p = as.numeric(aux[3,])
  
  nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados Covariaveis/resultado", i,"/hsd", sep="")
  
  aux = read.table(nome, sep = "%")
  aux = as.matrix(aux)
  
  aux = unlist(strsplit(aux, split = "Error:"))
  
  aux[1:10]
  
  aux = matrix(aux, nrow = 2)
  
  hsd = as.numeric(aux[2,])
  
  nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados Covariaveis/resultado", i,"/sexo", sep="")
  
  aux = read.table(nome, sep = "%")
  aux = as.matrix(aux)
  
  aux = unlist(strsplit(aux, split = "p ="))
  aux = unlist(strsplit(aux, split = "(Not", fixed = TRUE))
  aux = unlist(strsplit(aux, split = "(S", fixed = TRUE))
  
  aux[1:10]
  
  aux = matrix(aux, nrow = 3)
  
  psexo = as.numeric(aux[2,])
  
  nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados Covariaveis/resultado", i,"/idade", sep="")
  
  aux = read.table(nome, sep = "%")
  aux = as.matrix(aux)
  
  aux = unlist(strsplit(aux, split = "p ="))
  aux = unlist(strsplit(aux, split = "(Not", fixed = TRUE))
  aux = unlist(strsplit(aux, split = "(S", fixed = TRUE))
  
  aux[1:10]
  
  aux = matrix(aux, nrow = 3)
  
  pidade = as.numeric(aux[2,])
  
  class(cn)

  resultado = cbind(h,p,cn,pcn,pidade,psexo)
  resultado = resultado[order(resultado[,3]),]
  
  
 
  
  aux = read.table(nome)
  aux = aux[,]
  
  resultado = cbind(aux,resultado)
  
  head(resultado)
  dim(resultado)
  
  colnames(resultado) = c("Chr", "Start", "End", "H2r", "p-value", "CN", "pCN", "pAge", "pSex")
  
  final = rbind(final, resultado)
  print(i)
}

final = as.data.frame(final)
colnames(final) = c("Chr", "Start", "End", "H2r", "p", "CN", "pCN", "pAge", "pSex")
attach(final)

summary(final)
head(final)
dim(final)
tail(final)

#install.packages("qqman")
require(qqman)

pos = apply(final[,2:3], 1, mean)

nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados Covariaveis/man","pH",sep="")

jpeg(nome)
manH = as.data.frame(cbind(CN,Chr,pos,p))
colnames(manH)=c("SNP", "CHR", "BP","P")  
head(manH)
manhattan(manH, col = c("blue", "orange"))
dev.off()

nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados Covariaveis/man","pCN",sep="")

jpeg(nome)
manH = as.data.frame(cbind(CN,Chr,pos,pCN))
colnames(manH)=c("SNP", "CHR", "BP","P")  
head(manH)
manhattan(manH, col = c("blue", "orange"))
dev.off()


nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados Covariaveis/man","pAge",sep="")
jpeg(nome)
manH = as.data.frame(cbind(CN,Chr,pos, pAge))
colnames(manH)=c("SNP", "CHR", "BP","P")  
head(manH)
manhattan(manH, col = c("blue", "orange"))
dev.off()

nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados Covariaveis/man","pSex",sep="")
jpeg(nome)
manH = as.data.frame(cbind(CN,Chr,pos, pSexo))
colnames(manH)=c("SNP", "CHR", "BP","P")  
head(manH)
manhattan(manH, col = c("blue", "orange"))
dev.off()

summary(pCN)

length(which(pCN==0))

final[which(pCN==0),]

########## Selecionar resultados da regressao altura cn ##########

par(mfrow=c(1,1))

final = c()

for(i in 2:22){
  nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados CN/resultado", i,"/copy", sep="")
  
  aux = read.table(nome, sep = "%")
  aux = as.matrix(aux)
  
  aux = unlist(strsplit(aux, split = "cn."))
  aux = unlist(strsplit(aux, split = "p ="))
  aux = unlist(strsplit(aux, split = "(Not", fixed = TRUE))
  aux = unlist(strsplit(aux, split = "(S", fixed = TRUE))
  
  aux[1:10]
  
  aux = matrix(aux, nrow = 4)
  
  cn = as.numeric(aux[2,])
  cn[which(is.na(cn))] = 0 
  pcn = as.numeric(aux[3,])
  
  nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados CN/resultado", i,"/hp", sep="")
  
  aux = read.table(nome, sep = "%")
  aux = as.matrix(aux)
  
  aux = unlist(strsplit(aux, split = "H2r is"))
  aux = unlist(strsplit(aux, split = "p ="))
  aux = unlist(strsplit(aux, split = "(Not", fixed = TRUE))
  aux = unlist(strsplit(aux, split = "(S", fixed = TRUE))
  
  aux[1:10]
  
  aux = matrix(aux, nrow = 4)
  
  h = as.numeric(aux[2,])
  p = as.numeric(aux[3,])
  
  nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados CN/resultado", i,"/hsd", sep="")
  
  aux = read.table(nome, sep = "%")
  aux = as.matrix(aux)
  
  aux = unlist(strsplit(aux, split = "Error:"))
  
  aux[1:10]
  
  aux = matrix(aux, nrow = 2)
  hsd = as.numeric(aux[2,])
  
  resultado = cbind(h,p,cn,pcn)
  resultado = resultado[order(resultado[,3]),]
  
  
  nome = paste("/home/cicconella/Dropbox/2016/files/positions5", i, sep="")
  
  aux = read.table(nome)

  resultado = cbind(aux,resultado)
  
  head(resultado)
  dim(resultado)
  
  colnames(resultado) = c("Chr", "Start", "End", "H2r", "p-value", "CN", "pCN")
  
  final = rbind(final, resultado)
  print(i)
}

final = as.data.frame(final)
colnames(final) = c("Chr", "Start", "End", "H2r", "p", "CN", "pCN")
attach(final)

summary(final)
head(final)
dim(final)
tail(final)

#install.packages("qqman")
require(qqman)

pos = apply(final[,2:3], 1, mean)

nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados CN/man","pH",sep="")

jpeg(nome)
manH = as.data.frame(cbind(CN,Chr,pos,p))
colnames(manH)=c("SNP", "CHR", "BP","P")  
head(manH)
manhattan(manH, col = c("blue", "orange"))
dev.off()

nome = paste("/home/cicconella/Dropbox/2016/Altura/Resultados CN/man","pCN",sep="")

jpeg(nome)
manH = as.data.frame(cbind(CN,Chr,pos,pCN))
colnames(manH)=c("SNP", "CHR", "BP","P")  
head(manH)
manhattan(manH, col = c("blue", "orange"))
dev.off()

which(-log10(pCN)==max(-log10(pCN)))
final[which(-log10(pCN)==max(-log10(pCN))),]
