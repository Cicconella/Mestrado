
ana1 = read.table("/Users/Ana/Downloads/ana.ped", header = T, sep=",")
head(ana1)
dim(ana1)


ana2 = read.table("/Users/Ana/Downloads/ana.phen", header = T, sep=",")
head(ana2)
dim(ana2)



nubia1 = read.table("/Users/Ana/Downloads/ex.ped", header = T, sep=",")
head(nubia1)
dim(nubia1)


nubia2 = read.table("/Users/Ana/Downloads/ex.phen", header = T, sep=",")
head(nubia2)
dim(nubia2)


ped = cbind(info$IID,info$PAT,info$MAT,info$SEX,info$FID,info$Idade,info$altura)
colnames(ped) = c("id","fa","mo","sex","fid","idade","altura")

write.table(ped, "/home/cicconella/Dropbox/2016/ana2016-2.ped", row.names = F, quote = F, sep = ",")

########## Criando arquivo com os nomes dos fenotipos ##########

for(i in 1:22){
  nome = paste("/home/cicconella/Dropbox/2016/files/5phen",i,".phen", sep="")
  phen = read.table(nome, sep = ",", header = TRUE)
  aux = colnames(phen)[-c(1,2,3,4)]
  
  nome = paste("/home/cicconella/Dropbox/2016/files/5nomes",i, sep="")
  write.table(aux, nome, row.names = F, col.names = F, quote = F)
}
