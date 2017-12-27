require(kinship2)

##### Information about CNV regions #####
cnv = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/tableCNV", sep="\t", header = F)
colnames(cnv) = c("Chr","Start","End","Number","Length", "State", "CN", "Sample", "First Marker", "Last Marker")
head(cnv)
summary(cnv)
dim(cnv)
attach(cnv)

length(unique(Sample))

##### Information about quality control PennCNV per sample #####
qc = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/tableQC", sep="\t", header = T)
qc = qc[order(qc$File),] 
head(qc)

stargazer::stargazer(head(qc), rownames = F, summary = F)

summary(qc)
dim(qc)

########## Cleaning Bad Samples ##########


### LRR_SD

png("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/plots/lrr_sd.png")
hist(qc$LRR_SD, nc = 100, col="blue", xlab = "LRR Standard Deviation", ylim = c(0,70),
     main="Distribution of LRR Standard Deviation")
lines(rep(0.35,100),seq(0,70,by=(length(0:70)-1)/99))
dev.off()

length(which(qc$LRR_SD > 0.35))

### BAF_Mean


length(which(qc$BAF_mean > 0.6))+length(which(qc$BAF_mean < 0.4))

png("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/plots/baf_mean.png")
hist(qc$BAF_mean, col="blue", xlab = "BAF Mean", xlim =c(0.49,0.51),
     main = "Distribution of BAF Mean" )
dev.off()

### BAF_Drif

length(which(qc$BAF_drift > 0.01))

png("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/plots/baf_drift.png")
hist(qc$BAF_drift, nc = 100, col="blue", xlab = "BAF Drift", 
     main = "Distribution of BAF Drift")
lines(rep(0.01,100),seq(0,80,by=(length(0:80)-1)/99))
dev.off()

### WF

png("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/plots/wf.png")
hist(qc$WF, nc = 100, col="blue", xlab = "Waviness Factor", 
     main = "Distribution of Waviness Factor" )
lines(rep(0.04,100),seq(0,160,by=(length(0:160)-1)/99))
lines(rep(-0.04,100),seq(0,160,by=(length(0:160)-1)/99))
dev.off()

length(which(qc$WF > 0.04))+length(which(qc$WF < -0.04))

fail1 = which(qc$LRR_SD > 0.35)
fail1
fail2 = which(qc$BAF_drift > 0.01)
fail2
fail3 = which(qc$WF > 0.04)
fail3
fail4 = which(qc$WF < -0.04)
fail4


fail = unique(c(fail1, fail2, fail3, fail4))

length(fail)


qc = qc[-fail,]
dim(qc)
head(qc)
attach(qc)

cnv = merge(cnv, qc, by.x = "Sample", by.y = "File")
cnv = cnv[,-c(11:19)]
dim(cnv)
head(cnv)
summary(cnv)
attach(cnv)
ind = unique(Sample)
length(ind)

##### Individuos #####

infos = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/informacoesIndividuos", header=T)
head(infos)

#Classificados errados pelo PennCNV
errados = read.table("/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/errados")

dim(infos)
length(unique(infos$FID))

#Apenas os genotipados
a = infos[which(infos$GENOTIPADO==1),]
dim(a)
head(a)
length(unique(a[,1])) #familias

#Com qualidade
head(a)
b = merge(a, qc, by.x = "cel", by.y = "File")
head(b)
dim(b)
length(unique(b[,2]))

#So para checar quais são os individuos que sairam e tinha problema na estimativa do sexo
t = c()
n = c()
for(i in errados[,1]){
  bla = length(which(b$cel==i))
  if(bla==0){
    n = c(n,i)
  }else{
    t = c(t,i)
  }
}
t
n

head(b)
which(b$cel==1372)


c = merge(a, as.data.frame(n), by.x = "cel", by.y = "n")
dim(c)
c

#Retira familias sem relacao

dim(b)
head(b)
length(unique(b[,2]))

fam = unique(b[,2])

table(b$FID)

infos[infos$FID==119,]

head(infos)

i=87
semR =c()
for(i in fam){
  kin = infos[infos$FID==i,c(1,2,3,4,5,6)]
  
  head(kin)  
  
  kin$GENOTIPADO[which(is.na(kin$GENOTIPADO))]=0
  temG = which(kin$GENOTIPADO==1)
  
  kin <- with(kin, pedigree(kin$IID, kin$PAT, kin$MAT, kin$SEX, famid = kin$FID, affected = kin$altura))
  kin
  
  bla = round(8*kinship(kin))/4
  
  
  if(sum(apply(bla[temG,temG],1,sum)) == nrow(bla[temG,temG]))
    semR = c(semR,i)
}

semR


dim(b)
head(b)

d=b
for(i in semR)
  d = d[-which(d$FID==i),]


dim(d)
head(d)
tail(d)
length(unique(d[,2]))
sort(unique(d[,2]))

write.table(d, "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/Mestrado/ind_limpos",
            row.names = F)

head(d)
dim(d)

semR

infos[infos$FID==119,]


