
info = read.table("/home/cicconella/DadosMestrado/dados2", header = T, sep = ",")

head(info)  
summary(info)
dim(info)
str(info)
attach(info)


length(which(GENOTIPADO==1))


# Association between sample id and celfiles

ind = read.table("/home/cicconella/DadosMestrado/individuos")
colnames(ind) = c("cel", "IID")

class(ind)
head(ind)
dim(ind)
summary(ind)
str(ind)

# Associacao 

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


write.table(info, "/home/cicconella/DadosMestrado/informacoesIndividuos", row.names = F)



celfiles = cel[-which(is.na(cel))]

length(celfiles)

head(info)

summary(info)

length(which(GENOTIPADO==0))

a = info[which(GENOTIPADO==1),9]

length(a)

sort(a)[100:105]

write.table(sort(a), "/home/cicconella/DadosMestrado/celfiles", row.names = F)
