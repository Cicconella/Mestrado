##### Dados CNV #####

dir = "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado"

a = read.table(paste(dir,"/tableCNV", sep=""))
head(a)

colnames(a) = c("Chr", "Begin", "End", "NumberSNPs", "Size", "State", "CN",
                "Sample", "FirstMarker", "LastMarker")

##### Individuos Limpos #####

ind = read.table(paste(dir,"/ind_limpos", sep=""), header=T)
head(ind)

individuos = ind$cel
individuos = as.data.frame(individuos)
head(individuos)
head(a)
dim(a)

a = merge(a, individuos, by.x = "Sample", by.y = "individuos")
head(a)
dim(a)

attach(a)
head(a)


##### Samples #####
Sample
print(paste("Total de amostras:", length(unique(Sample))))

summary(as.numeric(table(Sample)))
sd(as.numeric(table(Sample)))

png(paste(dir, "/numberCNVs.png", sep=""))
hist(table(Sample), nc = 1000, main = "Number of CNVs per sample", xlab = "Number of CNVs")
dev.off()

table(Sample)[which(table(Sample)>500)]
paste(length(which(table(Sample)>500))/length(unique(Sample)), "são maiores que 500")
paste(length(which(table(Sample)>250))/length(unique(Sample)), "são maiores que 250")
paste(length(which(table(Sample)>100))/length(unique(Sample)), "são maiores que 100")

table(Sample)
summary(as.numeric(table(Sample)))

summary(table(Sample))

porcentagens = c()

for(i in seq(0,4000, by =25)){
  porcentagens = c(porcentagens, length(which(table(Sample)>i))/length(unique(Sample)))
}

porcentagens


png(paste(dir, "/samplesSize.png", sep=""))
plot(seq(0,4000, by =25), porcentagens, pch="", ylab="% Samples", xlab = "Number of CNVs")
lines(seq(0,4000, by =25), porcentagens, main = "Frequency of samples by number of CNVs")
dev.off()


cbind(seq(0,4000, by =25), porcentagens)

print(paste("% de pessoas com menos de 100 CNVs:", 1-porcentagens[5]))

##### Size

Size
summary(Size)

hist(Size, nc=100)

mean(Size)
medias = aggregate(Size, by = list(Sample), FUN=mean)

png(paste(dir, "/sizeCNVs.png", sep=""))
barplot(medias[,2], xlab = "Sample", ylab = "Mean length of CNVs")
dev.off()

as.numeric(which(table(Sample)>75))


png(paste(dir, "/lenghtByNumber.png", sep=""))
plot(as.numeric(table(Sample)),medias[,2], pch=".", xlab = "Number of CNVs",
     ylab = "Mean length of CNVs")
dev.off()


##### By chromosome

png(paste(dir, "/chrCNV.png", sep=""))
plot(table(Chr), xlab = "Chromosome", ylab = "Number of CNVs", cex.axis=0.75, 
     main = "CNV regions per Chromosome")
dev.off()

minimos = aggregate(Begin, by = list(Chr), FUN=min)
maximos = aggregate(Begin, by = list(Chr), FUN=max)

tamanhos = read.table(paste(dir,"/tamanhos", sep=""))

cbind(tamanhos, maximos)



