a = read.table("/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/Mestrado/CNV-Est/tableCNV")

head(a)

colnames(a) = c("Chr", "Begin", "End", "NumberSNPs", "Size", "State", "CN",
                "Sample", "FirstMarker", "LastMarker")
attach(a)
head(a)

##### Samples #####
Sample
print(paste("Total de amostras:", length(unique(Sample))))

summary(as.numeric(table(Sample)))
sd(as.numeric(table(Sample)))

png("/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/Mestrado/CNV-Est/numberCNVs.png")
hist(table(Sample), nc = 1000, main = "Number of CNVs per sample", xlab = "Number of CNVs")
dev.off()

table(Sample)[which(table(Sample)>500)]
paste(length(which(table(Sample)>500))/length(unique(Sample)), "são maiores que 500")
paste(length(which(table(Sample)>250))/length(unique(Sample)), "são maiores que 250")
paste(length(which(table(Sample)>100))/length(unique(Sample)), "são maiores que 100")

porcentagens = c()

for(i in seq(0,4000, by =25)){
  porcentagens = c(porcentagens, length(which(table(Sample)>i))/length(unique(Sample)))
}

porcentagens

png("/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/Mestrado/CNV-Est/samplesSize.png")
plot(seq(0,4000, by =25), porcentagens, pch="", ylab="% Samples", xlab = "Number of CNVs")
lines(seq(0,4000, by =25), porcentagens, main = "Frequency of samples by number of CNVs")
dev.off()

##### Size

Size
summary(Size)

hist(Size, nc=100)

mean(Size)
medias = aggregate(Size, by = list(Sample), FUN=mean)

png("/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/Mestrado/CNV-Est/sizeCNVs.png")
barplot(medias[,2], xlab = "Sample", ylab = "Mean length of CNVs")
dev.off()

as.numeric(which(table(Sample)>75))

png("/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/Mestrado/CNV-Est/lenghtByNumber.png")
plot(as.numeric(table(Sample)),medias[,2], pch=".", xlab = "Number of CNVs",
     ylab = "Mean length of CNVs")
dev.off()


##### By chromosome

png("/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/Mestrado/CNV-Est/chrCNV.png")
plot(table(Chr), xlab = "Chromosome", ylab = "Number of CNVs", cex.axis=0.75, 
     main = "CNV regions per Chromosome")
dev.off()

minimos = aggregate(Begin, by = list(Chr), FUN=min)
maximos = aggregate(Begin, by = list(Chr), FUN=max)

tamanhos = read.table("/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/Mestrado/CNV-Est/tamanhos")

cbind(tamanhos, maximos)

source("https://bioconductor.org/biocLite.R")
biocLite("Gviz")
