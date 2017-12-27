dir = "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado"

ind = read.table(paste(dir,"/ind_limpos", sep=""), header=T)
head(ind)
dim(ind)

table(ind$SEX)

length(unique(ind$FID))
summary(as.numeric(table(ind$FID)))

hist(ind$Idade, col="blue",nc=100)

table(ind$Idade,ind$SEX)

boxplot(Idade~SEX, data = rbind(transform(ind, SEX="0"), ind), col="blue", pch=16, main="Age Distribution", ylab="Age",
        names=c("All", "Male", "Female"))
summary(ind$Idade)


boxplot(altura~SEX, data = rbind(transform(ind, SEX="0"), ind), col="blue", pch=16, main="Height Distribution", ylab="Height",
        names=c("All", "Male", "Female"))
summary(ind$altura)
