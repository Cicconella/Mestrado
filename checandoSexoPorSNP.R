a = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado-Antigo/saida_snp", header = F)
head(a)


b = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado-Antigo/Antigo/saida2", header = F)
b[1:10, 1:10]

c = unlist(strsplit(as.character(b$V1), split ="-*-"))
c[1:30]
snps = c[seq(1,length(c), by=3)]
id = c[seq(2,length(c), by=3)]

snps[1:10]
id[1:10]

snps = cbind(snps,id)
head(snps)



b = cbind(apply(snps, 1, paste, collapse ="-"), c[seq(3,length(c), by=3)], b)
colnames(b)[1:2] = c("SNP", "Alelo")
b[1:10, 1:10]

head(a)
ref = a[,1:3]
colnames(ref) = c("SNP", "Chr", "Pos")
head(ref)

class(b)
class(ref)

bla = merge(ref, b, by.x="SNP", by.y = "SNP")
  bla[1:10,1:10]


table(bla$Chr)

?aggregate


xa = (intersect(which(bla$Chr=="X"),which(bla$Alelo=="A")))
xb = (intersect(which(bla$Chr=="X"),which(bla$Alelo=="B")))
ya = (intersect(which(bla$Chr=="Y"),which(bla$Alelo=="A")))
yb = (intersect(which(bla$Chr=="Y"),which(bla$Alelo=="B")))

bla[1:10,1:10]

teste = bla[,6]

ponto = c()


nova = bla[,-c(1:5)]
nova[1:10,1:10]

for(i in 1:nrow(nova)){
  nova[i,] = nova[i,]/max(nova[i,])  
}


for(i in 1:ncol(nova)){
  
  teste = nova[,i]
  
  um = (teste[c(xa)])
  dois = (teste[c(xb)])
  tres = (teste[c(ya)])
  quatro = (teste[c(yb)])
  
  um = mean(um)
  dois = mean(dois)
  tres = mean(tres)
  quatro = mean(quatro)
  
  ponto = rbind(ponto, c(mean(c(um, dois)), mean(c(tres,quatro))))
  ponto
  
}





sex = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado-Antigo/Antigo/file_sex.txt", 
                 header = T)

head(sex)

sex$computed_gender = as.character(sex$computed_gender)

azul = rgb(0,0,1,alpha=0.5)
vermelho = rgb(1,0,0,alpha=0.5)

sex[which(sex[,2]=="male"),2] = azul
sex[which(sex[,2]=="female"),2] = vermelho


plot(ponto[,1], ponto[,2], pch = 16, col = sex$computed_gender, xlab = "Mean X Chromosome Intensity (9508 SNPs)",
     ylab = "Mean Y Chromosome Intensity (134 SNPs)", xlim = c(0.7,0.9), ylim = c(0.7,0.9))
