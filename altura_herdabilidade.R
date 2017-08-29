library(kinship2)
library(coxme)

dados_limpos = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/ind_limpos", header=T)
head(dados_limpos)

dim(dados_limpos)

boxplot(dados_limpos$altura~dados_limpos$SEX, names = c("Men", "Women"))

table(dados_limpos$SEX)


head(dados_limpos)

dados = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/informacoesIndividuos", header=T)

head(dados)

table(dados$SEX)

pedig = with(dados, pedigree(id=dados$IID, dadid=dados$PAT, momid=dados$MAT, 
                           sex=(dados$SEX), famid=dados$FID, missid=0))
kmat = kinship(pedig)

fit <- lmekin(altura~1+Idade+SEX+(1|IID), data=dados, varlist=2*kmat,vinit=2)

varErro = fit$sigma^2                  
varPol = unlist(fit[3])
varTotal = varErro + varPol

hPol = varPol/varTotal
valor = hPol  

head(dados)

##### Com CNV #####

j=9
for(j in 16:22){
  
  cnv = read.table(paste("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/cromo", j, sep=""), header = T)
  cnv[1:10,1:10]
  
  x = unlist(strsplit(colnames(cnv)[-c(1:3)], split="X"))
  x = x[seq(2,length(x), by=2)]
  colnames(cnv)[-c(1:3)] = x
  
  head(dados)
  
  dim(cnv)
  
  h = rep(1,nrow(cnv))
  
  for(i in 1:length(h)){
    cn = cnv[i,-c(1:3)]
    
    cn = cbind(colnames(cn),as.numeric(as.character(cn)))
    cn = as.data.frame(cn)
    head(cn)
    summary(cn)
    
    dadoscn = merge(dados, cn, by.x="cel", by.y="V1", all.x = T)
    
    head(dadoscn)
    
    fit <- lmekin(altura~1+Idade+SEX+V2+(1|IID), data=dadoscn, varlist=2*kmat,vinit=2)
    
    varErro = fit$sigma^2                  
    varPol = unlist(fit[3])
    varTotal = varErro + varPol
    
    hPol = varPol/varTotal
    h[i] = hPol
    print(i)
  }
  
  png(paste("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/chr",j,".png", sep=""))
  plot(h, pch=16, ylim=c(0,1), main = paste("Chromosome", j))
  abline(valor,0, col = "red")
  dev.off()
    
}

