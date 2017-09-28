##### Criando .ped

dados = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/informacoesIndividuos", header=T)
head(dados)
dim(dados)

dim(dados)

ped = cbind(dados$IID,dados$PAT,dados$MAT,dados$SEX,dados$FID)
colnames(ped) = c("id","fa","mo","sex","fid")

head(ped)
dim(ped)

write.table(ped, "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/samples.ped", row.names = F, quote = F, sep = ",")
write.table(ped, "/home/cicconella/solar831/Chr1/samples.ped", row.names = F, quote = F, sep = ",")


##### Criando .phen para altura

head(dados)

phen = cbind(dados$IID, dados$altura, dados$Idade)
colnames(phen) = c("id","altura", "idade")
head(phen)
dim(phen)

phen[which(is.na(phen))] = ""

write.table(phen, "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/samples.phen", row.names = F, quote = F, sep = ",")
write.table(phen, "/home/cicconella/solar831/Chr1/samples.phen", row.names = F, quote = F, sep = ",")

##### Criando .phen para CNV

j=2
for(j in 4:22){
  cnv = read.table(paste("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/cromo", j, sep=""), header = T)
  dim(cnv)
  
  cnv[1:10,1:10]
  
  x = unlist(strsplit(colnames(cnv)[-c(1:3)], split="X"))
  x = x[seq(2,length(x), by=2)]
  colnames(cnv)[-c(1:3)] = x
  
  dim(cnv)
  
  cnv[1:10,1:10]
  
  
  cn = cnv[,-c(1:3)]
  cn[1:10,1:10]
  cn = cbind(colnames(cn),t(as.matrix(cn)))
  cn[1:10,1:10]
  cn = as.data.frame(cn)
  cn[1:10,1:10]
  
  dadoscn = merge(dados, cn, by.x="cel", by.y="V1", all.x = T)
  dadoscn[1:10,1:10]
  dadoscn = dadoscn[,-c(1,2,4,5,6,7)]
  
  dadoscn[1:10,1:10]
  
  for(i in 4:ncol(dadoscn)){
    colnames(dadoscn)[i] = paste("CN",(i-3), sep="")  
  }
  
  colnames(dadoscn)[1]= "id"
  
  dadoscn[1:10,1:10]

  colnames(dadoscn)
    
  
  write.table(dadoscn, paste("/home/cicconella/solar831/Chr",j,"/cn.phen", sep=""), row.names = F, quote = F, sep = ",", na="")
  write.table(colnames(dadoscn)[-c(1:3)], paste("/home/cicconella/solar831/Chr",j,"/files", sep=""), row.names = F, quote = F, sep = ",", na="", col.names = F)

}

##### Criando .phen para CNV

j=1

for(j in 1:22){
  a = read.table(paste("/home/cicconella/solar831/Chr",j,"/cn.phen", sep=""), header = T, sep=",")
  class(a)
  
  a[1:10,1:10]
  
  bla = a[,-c(1:3)]
  
  bla[1:10,1:10]
  
  bla[bla==2] = "Normal"
  bla[1:10,1:10]
  bla[bla!="Normal"] = 0
  bla[1:10,1:10]
  bla[bla=="Normal"] = 1
  bla[1:20,1:20]
  
  a[,-c(1:3)] = bla
  a[1:10,1:10]
  
  write.table(a, paste("/home/cicconella/solar831/Chr",j,"/cnBin.phen", sep=""), row.names = F, quote = F, sep = ",", na="")
  
}
head(a)
