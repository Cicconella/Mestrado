a = "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/Cromossomos/cromo9"

a = read.table(a, header = T)

dim(a)
a[1:10,1:10]

x = a[which(a[,2]==78960219):which(a[,3]==78967224),]
x[,1:10]

table(as.numeric(x[1,-c(1:3)]))
table(as.numeric(x[2,-c(1:3)]))

dados_limpos = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/ind_limpos", header=T)
head(dados_limpos)

x = t(x)
head(x)
tail(x)

which(x[,1]==3)

n = unlist(strsplit(rownames(x), split="X"))
n[1:20]
n=n[c(1:3,seq(5,length(n),by=2))]
dim(x)
length(n)
x = cbind(x,n)
head(x)

x = as.data.frame(x)
head(x)
table(x$`423`)

head(dados_limpos)

fim = merge(x,dados_limpos[,c(1,6,8,9)], by.x="n", by.y="cel")
head(fim)

table(fim[which(fim$`423`==3),"SEX"])

fim$`423` = as.character(fim$`423`)
fim$`423` = as.factor(fim$`423`)

fim$`424` = as.character(fim$`424`)
fim$`424` = as.factor(fim$`424`)

boxplot(fim$altura~fim$`423`, xlab = "Number of Copies", ylab="Height", col="blue")
boxplot(fim$altura~fim$`424`, xlab = "Number of Copies", ylab="Height", col="blue")
