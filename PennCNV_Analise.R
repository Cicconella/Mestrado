
##### Comparando Log R Ratio #####

novo = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/rb.0002", header = T, sep="\t")
antigo = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado-Antigo/rb.1186ACP0002", header = T, sep="\t")

colnames(novo) = c("Name", "Chr", "Pos", "LRR", "BAF")
colnames(antigo) = c("Name", "Chr", "Pos", "LRR", "BAF")

head(novo)
head(antigo)

dim(novo)
dim(antigo)

dim(novo)[1]-dim(antigo)[1]

ambos = merge(novo, antigo, by.x = "Name", by.y = "Name", all.x = T)

head(ambos)
summary(ambos)
attach(ambos)

length(which(Chr.x==Chr.y))
length(which(Pos.x==Pos.y))
length(which(LRR.x==LRR.y))
length(which(BAF.x==BAF.y))


barplot(table(Chr.x))

##### Comparando CNV #####

cnv = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/tableCNV", header = F, sep="\t")
head(cnv)

