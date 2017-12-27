library(kinship2)
#citation("kinship2")
library(coxme)
library(qqman)

extract_coxme_table <- function (mod){
  beta <- mod$coefficients$fixed
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}


dados_limpos = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/ind_limpos", header=T)
head(dados_limpos)

dim(dados_limpos)

boxplot(dados_limpos$altura~dados_limpos$SEX, names = c("Men", "Women"))

table(dados_limpos$SEX)


head(dados_limpos)

dados = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/informacoesIndividuos", header=T)

head(dados)
dim(dados)

ancestralidade = read.table("pcs-mariza-abril2017.csv", header = T, sep = ",")
head(ancestralidade)
dim(ancestralidade)

ancestralidade = ancestralidade[,c(1,6,7)]
head(ancestralidade)

dados = merge(dados, ancestralidade, by="IID", all.x=T)
head(dados)
dim(dados)

cor = rgb(0,0,1,0.5)

plot(ancestralidade[,c(2,3)], pch=20, main="Individual Ancestry", xlab="PC1",
     ylab="PC2", col=cor)

table(dados$SEX)

pedig = with(dados, pedigree(id=dados$IID, dadid=dados$PAT, momid=dados$MAT, 
                             sex=(dados$SEX), famid=dados$FID, missid=0))
kmat = kinship(pedig)

fit <- lmekin(altura~1+Idade+SEX+pca1+pca2+(1|IID), data=dados, varlist=2*kmat,vinit=2)

varErro = fit$sigma^2                  
varPol = unlist(fit[3])
varTotal = varErro + varPol

hPol = varPol/varTotal
valor = hPol  

head(dados)

valor


##### Com CNV #####

herdabilidades_dico = c()
herdabilidades_continuo = c()
herdabilidades_factor = c()

pv_dico = c()
pv_continuo = c()
pv_factor = list()

bcnv_dico = c()
bcnv_continuo = c()
bcnv_factor = list()

i=1
j=1

for(j in 1:22){
  
  cnv = read.table(paste("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/cromo", j, sep=""), header = T)
  cnv[1:10,1:10]
  dim(cnv)
  
  x = unlist(strsplit(colnames(cnv)[-c(1:3)], split="X"))
  x = x[seq(2,length(x), by=2)]
  colnames(cnv)[-c(1:3)] = x
  
  head(dados)
  
  dim(cnv)
  
  h_dico = rep(NA,nrow(cnv))
  h_continuo = rep(NA,nrow(cnv))
  h_fator = rep(NA,nrow(cnv))
  
  p_dico = rep(NA,nrow(cnv))
  p_continuo = rep(NA,nrow(cnv))
  p_fator = list()
  
  b_dico = rep(NA,nrow(cnv))
  b_continuo = rep(NA,nrow(cnv))
  b_fator = list()
  
  #nrow(cnv)
  for(i in 1:nrow(cnv)){
    
    cn = cnv[i,-c(1:3)]
    
    cn = cbind(colnames(cn),as.numeric(as.character(cn)))
    cn = as.data.frame(cn)
    head(cn)
    summary(cn)
    
    class(cn[,2])
    
    dadoscn = merge(dados, cn, by.x="cel", by.y="V1", all.x = T)
    
    head(dadoscn)
    
    class(cn$V2)
    table(cn$V2)
    
    head(dadoscn)
    table(dadoscn$V2)
    
    dadoscn <- within(dadoscn, V2 <- relevel(V2, ref = 2))
    
    table(dadoscn$V2)
    
    colnames(dadoscn)[colnames(dadoscn)=="V2"] = "CNV"
    
    dadoscn$CNV
    
    tryCatch({
      fit <- lmekin(altura~1+Idade+SEX+pca1+pca2+as.numeric(as.character(CNV))+(1|IID), data=dadoscn, varlist=2*kmat,vinit=2)
      #print(fit)
      #dadoscn$CNV[1:10]      b_fator[[i]] = b2

      #as.numeric(as.character(dadoscn$CNV[1:10]))
      fit
      class(fit)
      fit[1]
      fit$var
      
      p = extract_coxme_table(fit)
      p = p[6,4]
      
      p_continuo[i] = p
      
      varErro = fit$sigma^2                  
      varPol = unlist(fit[3])
      varTotal = varErro + varPol
      
      hPol = varPol/varTotal
      hPol
      h_continuo[i] = hPol
      
      b = extract_coxme_table(fit)
      b = b[6,1]

      b_continuo[i] = b
      
      }, error=function(e){
      h_continuo[i] = NA
      p_continuo[i] = NA
      b_continuo[i] = NA
      #print(i)
    }
    )
    
    tryCatch({
      fit <- lmekin(altura~1+Idade+SEX+pca1+pca2+CNV+(1|IID), data=dadoscn, varlist=2*kmat,vinit=2)
      #print(fit)
      #dadoscn$CNV[1:10]
      #as.numeric(as.character(dadoscn$CNV[1:10]))
      
      class(fit)
      fit
      
      p = extract_coxme_table(fit)
      p = p[6:nrow(p),4]
      
      p_fator[[i]] = p
      
      varErro = fit$sigma^2                  
      varPol = unlist(fit[3])
      varTotal = varErro + varPol
      
      hPol = varPol/varTotal
      hPol
      h_fator[i] = hPol
      #print(i)
      
      b = extract_coxme_table(fit)
      b2 = b[6:nrow(b),1]
      names(b2) = rownames(b)[6:nrow(b)] 
      
      b_fator[[i]] = b2
      
    }, error=function(e){
      h_fator[i] = NA
      p_fator[[i]] = NA
      b_fator[[i]] = NA
      
      #print(i)
    }
    )
    
    head(dadoscn)
    
    dadoscn$CNV = as.character(dadoscn$CNV)
    dadoscn$CNV[which(dadoscn$CNV==2)] = "X"
    dadoscn$CNV[which(dadoscn$CNV!="X")] = 1
    dadoscn$CNV[which(dadoscn$CNV=="X")] = 0
    
    table(dadoscn$CNV)
    class(dadoscn$CNV)
    
    dadoscn$CNV = as.factor(dadoscn$CNV)
    
    tryCatch({
      fit <- lmekin(altura~1+Idade+SEX+pca1+pca2+CNV+(1|IID), data=dadoscn, varlist=2*kmat,vinit=2)
      #print(fit)
      #dadoscn$CNV[1:10]
      #as.numeric(as.character(dadoscn$CNV[1:10]))
      
      class(fit)
      fit
      fit$var
      
      p = extract_coxme_table(fit)
      p = p[6,4]
      
      p_dico[i] = p
      
      varErro = fit$sigma^2                  
      varPol = unlist(fit[3])
      varTotal = varErro + varPol
      
      hPol = varPol/varTotal
      hPol
      h_dico[i] = hPol
      #print(i)
      
      b = extract_coxme_table(fit)
      b = b[6,1]
      
      b_dico[i] = b
      
    }, error=function(e){
      h_dico[i] = NA
      b_dico[i] = NA
      p_dico[i] = NA
      #print(i)
    }
    )
    
  }
  
  # png(paste("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/chr",j,".png", sep=""))
  # plot(h, pch=16, ylim=c(0,1), main = paste("Chromosome", j))
  # abline(valor,0, col = "red")
  # dev.off()
  
  h = cbind(rep(j, length(h_dico)), h_dico)
  herdabilidades_dico = rbind(herdabilidades_dico, h)
  head(herdabilidades_dico)
  
  h = cbind(rep(j, length(h_continuo)), h_continuo)
  herdabilidades_continuo = rbind(herdabilidades_continuo, h)
  head(herdabilidades_continuo)
  
  h = cbind(rep(j, length(h_fator)), h_fator)
  herdabilidades_factor = rbind(herdabilidades_factor, h)
  head(herdabilidades_factor)
  
  p = cbind(rep(j, length(p_dico)), p_dico)
  pv_dico = rbind(pv_dico, p)
  head(pv_dico)
  
  p = cbind(rep(j, length(p_continuo)), p_continuo)
  pv_continuo = rbind(pv_continuo, p)
  head(pv_continuo)
  
  pv_factor[[j]] = p_fator
  pv_factor[[j]]
  
  
  b = cbind(rep(j, length(b_dico)), b_dico)
  bcnv_dico = rbind(bcnv_dico, b)
  head(bcnv_dico)
  
  b = cbind(rep(j, length(b_continuo)), b_continuo)
  bcnv_continuo = rbind(bcnv_continuo, b)
  head(bcnv_continuo)
  
  bcnv_factor[[j]] = b_fator
  bcnv_factor[[j]]
  
  print(j)  
}



head(herdabilidades_dico)
head(herdabilidades_continuo)
head(herdabilidades_factor)
head(pv_dico)
summary(pv_dico)
which(pv_dico[,2]<0.001)
which(pv_dico[,2]<0.01)
which(pv_dico[,2]<0.05)
head(pv_continuo)
which(pv_continuo[,2]<0.001)
which(pv_continuo[,2]<0.01)
which(pv_continuo[,2]<0.05)

pv_factor



pv_fator = c()

for(j in 1:22){
  for(i in pv_factor[[j]]){
    pv_fator = c(pv_fator, min(i))
  }
}

bcnv_fator = c()
nomes=c()
j=1
i=bcnv_factor[[j]][[1]]

for(j in 1:22){
  for(i in bcnv_factor[[j]]){
    bcnv_fator = c(bcnv_fator, min(i))
    x = names(i[i==min(i)])
    if(length(x)>0){
      nomes = c(nomes,x)  
    }else{
      nomes = c(nomes,NA)
    }
     
  }
}

##### Resultados #####

info = c()

for(j in 1:22){
  
  cnv = read.table(paste("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/cromo", j, sep=""), header = T)
  cnv = cnv[,1:3]
  dim(cnv)
  head(cnv)
  
  info = rbind(info, cnv)
  
}

dim(info)
dim(pv_dico)


man = cbind(info, pv_dico)
man = as.data.frame(man)

head(man)

minimos = man
head(minimos)
dim(minimos)
dim(bcnv_dico)
head(bcnv_dico)
minimos = cbind(minimos, bcnv_dico)
head(minimos)
minimos = minimos[order(minimos$p_dico),]
minimos[1:10,]
minimos = minimos[,c(1,2,3,7,5)]
colnames(minimos) = c("Chromosome", "Start", "End", "Coefficient", "P-Value")
minimos[1:10,]
stargazer::stargazer(minimos[1:20,], summary = F, rownames = F)


table(man$chrom==man$V1)

man = cbind(man$chrom, apply(man[,2:3],1,mean),man$p_dico)
man = as.data.frame(man)
colnames(man) = c("CHR","BP", "P")
head(man)

png("man_height_dico_pcs.png", width = 1024, height = 452)
par(bg=NA)
manhattan(man, col = c("blue","black"), cex.main=2,
          main = "Height (CNV as Dichotomous Covariate)",
          cex.axis=0.75, ylim=c(0,5), suggestiveline = F)
dev.off()

dim(info)
dim(pv_continuo)

man = cbind(info, pv_continuo)
man = as.data.frame(man)

head(man)
table(man$chrom==man$V1)

minimos = man
head(minimos)
dim(bcnv_continuo)
head(bcnv_continuo)
minimos = cbind(minimos, bcnv_continuo)
head(minimos)
minimos = minimos[order(minimos$p_continuo),]
minimos[1:10,]
minimos = minimos[,c(1,2,3,7,5)]
colnames(minimos) = c("Chromosome", "Start", "End", "Coefficient", "P-Value")
minimos[1:10,]
minimos$`P-Value` = as.character(minimos$`P-Value`)
stargazer::stargazer(minimos[1:20,], summary = F, rownames = F)


man = cbind(man$chrom, apply(man[,2:3],1,mean),man$p_continuo)
man = as.data.frame(man)
colnames(man) = c("CHR","BP", "P")
head(man)

png("man_height_continous_pcs.png", width = 1024, height = 452)
par(bg=NA)
manhattan(man, col = c("blue","black"), cex.main=2,
          main = "Height (CNV as Continous Covariate)",
          cex.axis=0.75, ylim= c(0,5),suggestiveline = F)
dev.off()

dim(info)
length(pv_fator)
length(nomes)

man = cbind(info, pv_fator)
man = as.data.frame(man)

head(man)
#table(man$chrom==man$V1)


minimos = man
head(minimos)
minimos = cbind(minimos, nomes, bcnv_fator)
head(minimos)
minimos = minimos[order(minimos$pv_fator),]
minimos[1:10,]
minimos = minimos[,c(1,2,3,5,6,4)]
minimos[1:10,]
colnames(minimos) = c("Chromosome", "Start", "End","Factor","Coefficient", "P-Value")
minimos[1:10,]
minimos$`P-Value` = as.character(minimos$`P-Value`)
stargazer::stargazer(minimos[1:20,], summary = F, rownames = F)



man = cbind(man$chrom, apply(man[,2:3],1,mean),man$pv_fator)
man = as.data.frame(man)
colnames(man) = c("CHR","BP", "P")
head(man)

png("man_height_cat_pcs.png", width = 1024, height = 452)
par(bg=NA)
manhattan(man, col = c("blue","black"), cex.main=2,
          main = "Height (CNV as Categorical Covariate)",
          cex.axis=0.75, ylim=c(0,5), suggestiveline = F)
dev.off()

dim(info)
dim(herdabilidades_dico)
x = herdabilidades_dico

man = cbind(info, x)
man = as.data.frame(man)

head(man)
table(man$chrom==man$V1)

man[man$h_dico<valor,]

man = cbind(man$chrom, apply(man[,2:3],1,mean),man$h_dico)
man = as.data.frame(man)
colnames(man) = c("CHR","BP", "P")
head(man)

man[(which(man$P==min(man$P))-10):(which(man$P==min(man$P))+10),]

png("man_height_her_dico_pcs.png", width = 1024, height = 452)
par(bg=NA)
manhattan(man, col = c("blue","black"), cex.main=2,
          main = "Height Heritability (CNV as Dichotomous Covariate)",
          cex.axis=0.75, logp = F, ylim =c(0.8,0.9), suggestiveline = valor)
dev.off()

dim(info)
dim(herdabilidades_continuo)
x = herdabilidades_continuo

man = cbind(info, x)
man = as.data.frame(man)

head(man)
table(man$chrom==man$V1)

man[man$h_continuo<valor,]

man = cbind(man$chrom, apply(man[,2:3],1,mean),man$h_continuo)
man = as.data.frame(man)
colnames(man) = c("CHR","BP", "P")
head(man)

png("man_height_her_continous_pcs.png", width = 1024, height = 452)
par(bg=NA)
manhattan(man, col = c("blue","black"), cex.main=2,
          main = "Height Heritability (CNV as Continous Covariate)",
          cex.axis=0.75, logp = F, ylim =c(0.8,0.9), suggestiveline = valor)
dev.off()


dim(info)
dim(herdabilidades_factor)
x = herdabilidades_factor

man = cbind(info, x)
man = as.data.frame(man)

man[man$h_fator<valor,]

head(man)
table(man$chrom==man$V1)

man = cbind(man$chrom, apply(man[,2:3],1,mean),man$h_fator)
man = as.data.frame(man)
colnames(man) = c("CHR","BP", "P")
head(man)

png("man_height_her_cat_pcs.png", width = 1024, height = 452)
par(bg=NA)
manhattan(man, col = c("blue","black"), cex.main=2,
          main = "Height Heritability (CNV as Categorical Covariate)",
          cex.axis=0.75, logp = F, ylim =c(0.8,0.9), suggestiveline = valor)
dev.off()


