mil = read.table("dados1000", sep="\t")
dim(mil)
head(mil)
mil = mil[,-c(2,6:8)]
head(mil)
colnames(mil) = c("Chr", "Mut", "Start", "End", "Things")
head(mil)
attach(mil)

unique(Mut)
mil = mil[-which(mil$Mut=="alu_insertion"),]
mil = mil[-which(mil$Mut=="mobile_element_insertion"),]
head(mil)
dim(mil)
attach(mil)

#Lenght
summary(End-Start)

hist(log(End-Start), nc=100)

Things[1]

info = c()

for(i in 1:nrow(mil)){
  info = rbind(info, unlist(strsplit(as.character(Things[i]), split=";")))
}
i
dim(info)
info = as.data.frame(info)
head(info)

table(info$V1)
