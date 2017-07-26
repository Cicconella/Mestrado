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
summary(qc)
dim(qc)

########## Cleaning Bad Samples ##########
fail1 = intersect(which(qc$LRR_SD > 0.35), which(qc$BAF_drift > 0.01))
fail2 = intersect(which(qc$LRR_SD > 0.35), which(qc$WF > 0.04))
fail3 = intersect(which(qc$LRR_SD > 0.35), which(qc$WF < 0.04))

fail = c(fail2,fail3)

qc = qc[-fail,]
dim(qc)
head(qc)
attach(qc)
#  a = intersect(qc$File, cel_errados)

cnv = merge(cnv, qc, by.x = "Sample", by.y = "File")
cnv = cnv[,-c(11:19)]
dim(cnv)
head(cnv)
summary(cnv)
attach(cnv)
ind = unique(Sample)
length(ind)

remove(qc)
remove(fail3)
remove(fail2)
remove(fail1)
#remove(ind)


kjhiuhiuhiuh