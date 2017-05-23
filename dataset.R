########## Describing Dataset ##########

dataset = read.table("/home/cicconella/Dropbox/2016/Project/dataset", header = T)

head(dataset)
summary(dataset)

dataset2 = dataset[-which(dataset$GENOTIPADO==0),]
dataset2 = dataset2[-which(is.na(dataset2$GENOTIPADO)),]
head(dataset2)

# Numero de familias 

n = length(unique(dataset$FID))
n2 = length(unique(dataset2$FID))

