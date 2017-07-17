a = read.table("/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/Mestrado/Baependi.csv",
               header = T, sep=";")
colnames(a)
a = a[,c(1,2,3,4,5,195,196,197)]
head(a)
dim(a)


head(info)
dim(info)

a = merge(info, a, by.x = "IID", by.y = "ID")

head(a)

a[1:10,]

a$ALTURA[which(a$ALTURA=="*")] = NA
a$PESO[which(a$PESO=="*")] = NA
a$CIRCABD[which(a$CIRCABD=="*")] = NA

a = a[-which(a$ALTURA!=a$altura),]

head(a)

a = a[,c(1,2,3,4,5,6,7,8,14,16,9)]
head(a)

info = a

chromosomes[[1]][,1:10]
mutacoes = apply(chromosomes[[1]], 1, contaMut)
mutacoes[which(mutacoes==max(mutacoes))]

c = which(mutacoes==max(mutacoes))[2] 
genotipo = as.numeric(chromosomes[[1]][c,-c(1:3)])

table(genotipo)

genotipo = (cbind(colnames(chromosomes[[1]])[-c(1:3)], genotipo))

colnames(genotipo) = c("celfiles", "cnv")
genotipo = as.data.frame(genotipo)
head(genotipo)

head(info)

final = merge(info, genotipo, by.x = "cel", by.y = "celfiles")
final[1:100,]
dim(final)
attach(final)

final = final[-which(is.na(final$PESO)),]
final = final[-which(is.na(final$CIRCABD)),]

head(final)
final$PESO = as.numeric(as.character(final$PESO))
final$PESO = final$PESO/1000

final$PESO = as.numeric(as.character(final$PESO))
final$PESO = final$PESO/1000

final$altura = as.numeric(as.character(final$altura))
final$altura = final$altura/100
final$CIRCABD = as.numeric(as.character(final$CIRCABD))
final$cnv = as.numeric(as.character(final$cnv))

head(final)

imc = function(medidas){
  return(medidas[1]/(medidas[2]^2))
}

imcs = apply(final[,c(10,9)], 1, imc)

bla = cbind(imcs, final$cnv)

boxplot(imcs~final$cnv)
boxplot(final$CIRCABD~final$cnv)

table(final$SEX,final$cnv)










chromosomes[[1]][,1:10]
