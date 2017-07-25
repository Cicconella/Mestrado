##### Sexo baseado no celfile #####
sexo = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/file_sex_cel", sep="\t", header = T)
head(sexo)
cel = sexo$cel_files
cel = as.character(cel)
cel = substr(cel, 1,4)
sexo$cel_files = cel
colnames(sexo) = c("cel","sexo")
head(sexo)
attach(sexo)

##### Dados individuos #####
info = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/dados", header = T, sep = ",")
head(info)  
summary(info)
dim(info)
str(info)
attach(info)

##### CEL + IID #####
ind = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/cel+iid", header=T)
colnames(ind) = c("cel", "IID")
class(ind)
head(ind)
dim(ind)
summary(ind)
str(ind)


##### Associação dos arquivos #####
cel = IID 

for(i in 1:length(cel)){
  
  aux = which(ind$IID==cel[i])
  if(length(aux)!=0)
    cel[i] = max(ind$cel[aux])
  else
    cel[i] = NA
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

write.table(info, "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/informacoesIndividuos", row.names = F)

##### Conferindo com o identificado #####

head(info)
head(sexo)

sexo$cel = as.numeric(as.character(sexo$cel))

completo = merge(sexo, info, by.x = "cel", by.y = "cel")

head(completo)

table(completo[,c(2,7)])

cel_errados = c(completo[which(completo[,2]=="female" & completo[,7]==1),1],completo[which(completo[,2]=="male" & completo[,7]==2),1])


