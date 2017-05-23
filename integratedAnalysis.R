##### clear all variables ########
rm(list = ls())

##### Libraries #####
source("http://bioconductor.org/biocLite.R")
require(igraph)
require(plotrix)
library(plotly)
library(biomaRt)
library(STRINGdb)

##### Load expression information ######

### Insert the path to the folder (it will be use for loading and saving all the files)
folder = "/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/metaanalysis/Integrated Analysis/" 

### Load protein expression data
expP = read.table(paste(folder, "01d_PEAKS_Ensembl_1D25_unique.txt", sep = ""), header = T)
expP = cbind(expP[,1:2],expP[,4:11], expP$COPD10, expP$Control1, expP[,14:21], expP$Control10)

### Load diferential expression genes from proteomics analysis
proteins = read.table(paste(folder, "General_DEs_1D25.txt", sep = ""), header = T, sep =",")

### Combine information 
DEP = merge(proteins, expP, by.x = "Gene", by.y = "Gene")

rm(proteins, expP)

### Load transcripts expression data
expT = read.table(paste(folder ,"20160302_LundData_transcriptomics_freeze_exons.txt", sep = ""), header = T, sep = "\t")
attach(expT)

expT = as.data.frame(cbind(as.character(Gene), COPD1, COPD2, COPD3, COPD4, COPD5, COPD6, COPD7, COPD8, COPD9, COPD10,
                            Control1, Control2, Control3, Control4, Control6, Control7, Control8, Control10))

detach(expT)

for(i in 2:19){
  expT[,i] = as.numeric(as.character(expT[,i])) 
}

### Normalization: Upper Quantile or Counts per Million

### Upper quantile normalization
for(i in 2:ncol(expT)){
  a = as.numeric(as.character(expT[,i]))
  a = a/(quantile(a)[4])
  expT[,i] = a
}

### Counts per Million
# for(i in 2:ncol(expT)){
#   a = as.numeric(as.character(expT[,i]))
#   print(sum(a))
#   print(a[1:10])
#   b = sum(a)/1000000
#   a = a/b
#   print(a[1:10])
#   expT[,i] = a
# }

### Load diferential expression genes from transcriptomics analysis
proteins = read.table(paste(folder, "General_DEs_20160302_LundData_transcriptomics_freeze_exons.txt", sep = ""), header = T, sep =",")

DET = merge(proteins, expT, by.x = "X", by.y = "V1")

rm(proteins, expT, a, i)

head(DEP)
head(DET)

##### Networks - Separated #####

### Names of genes are in:
DEP$gene = as.character(DEP$gene)
DET$Associated.Gene.Name = as.character(DET$Associated.Gene.Name)
colnames(DET)[7] = "gene"

### Ploting networks based on String Database 

string_db = STRINGdb$new(version="10", species=9606, score_threshold=0.4, input_directory="")

genesP = string_db$map(DEP, "gene", removeUnmappedRows = TRUE)

rm(DEP)

#png(paste(folder,"stringNetworkP.png", sep = ""),  width = 1000, height = 1000)
#string_db$plot_network(genesP[,29])
#dev.off()

genesT = string_db$map(DET, "gene", removeUnmappedRows = TRUE)

rm(DET)

#png(paste(folder,"stringNetworkT.png", sep = ""),  width = 1000, height = 1000)
#string_db$plot_network(genesT[,27])
#dev.off()

### Selecting network with only active interaction sources from experimentals and calculate rank 

### For proteomics
links = string_db$get_interactions(genesP[,29])

links = links[links$experiments!=0,c(1,2)]
head(links)

# Change the String ID for gene names
for(i in 1:nrow(links)){
  links[i,1] = genesP[which(genesP$STRING_id==links[i,1]),1]
  links[i,2] = genesP[which(genesP$STRING_id==links[i,2]),1]
}

g = graph_from_data_frame(links, directed=F)

#png(paste(folder,"simpleNetworkP.png", sep = ""),  width = 1000, height = 1000)
plot(g, main = "Proteomics data", vertex.color = "green", vertex.label.dist=0.3, vertex.label.degree=pi/2, vertex.size=5, vertex.label.color = "gray")
#dev.off()

pr = page_rank(g)$vector
pr = cbind(names(pr),pr)
colnames(pr) = c("Gene", "Rank")
pr = as.data.frame(pr)
pr$Rank = as.numeric(as.character(pr$Rank))

head(genesP)

genesP = merge(genesP, pr, by.x = "gene", by.y = "Gene", all.x = T)
genesP$Rank[is.na(genesP$Rank)] = 0

### For transcripts
links = string_db$get_interactions(genesT[,27])

links = links[links$experiments!=0,c(1,2)]
head(links)

# Change the String ID for gene names
for(i in 1:nrow(links)){
  links[i,1] = genesT[which(genesT$STRING_id==links[i,1]),1]
  links[i,2] = genesT[which(genesT$STRING_id==links[i,2]),1]
}

g = graph_from_data_frame(links, directed=F)

#png(paste(folder,"simpleNetworkT.png", sep = ""),  width = 1000, height = 1000)
plot(g, main = "Transcripts data", vertex.color = "green", vertex.label.dist=0.3, vertex.label.degree=pi/2, vertex.size=5, vertex.label.color = "gray")
#dev.off()

pr = page_rank(g)$vector
pr = cbind(names(pr),pr)
colnames(pr) = c("Gene", "Rank")
pr = as.data.frame(pr)
pr$Rank = as.numeric(as.character(pr$Rank))

head(genesT)

genesT = merge(genesT, pr, by.x = "gene", by.y = "Gene", all.x = T)
genesT$Rank[is.na(genesT$Rank)] = 0


rm(links, pr, g, i)
rm(string_db)

##### Calculate weight for each gene #####

head(genesP)
head(genesT)

zeroP = 0 - min(genesP$logFC)
genesP$logFC = genesP$logFC-min(genesP$logFC)
zeroP = zeroP/sum(genesP$logFC) 
genesP$logFC = genesP$logFC/sum(genesP$logFC)

zeroT = 0 - min(genesT$logFC)
genesT$logFC = genesT$logFC-min(genesT$logFC)
zeroT = zeroT/sum(genesT$logFC) 
genesT$logFC = genesT$logFC/sum(genesT$logFC)

beta = 0.5

w = beta*genesP$Rank+(1-beta)*genesP$logFC
genesP = cbind(genesP, w)
head(genesP)

w = beta*genesT$Rank+(1-beta)*genesT$logFC
genesT = cbind(genesT, w)
head(genesT)

rm(w, beta)

##### Distances between samples #####

distance = function(a, b, w){
  
  c = sqrt(sum(((a-b)/(a+b+1))^2*w))  
  return(c)

}

aux  = dim(genesP)[2]-11

distances = matrix(NA, nrow = aux, ncol = aux)

for(i in 1:aux){
  for(j in i:aux){
    distances[i,j] = distance(genesP[,i+8],genesP[,j+8], genesP[,31])   
  }
}

#plot_ly(
#  x = seq(1:aux), y = seq(1:aux),
#  z = distances, type = "heatmap"
#)

aux  = dim(genesT)[2]-11

distances = matrix(NA, nrow = aux, ncol = aux)

for(i in 1:aux){
  for(j in i:aux){
    distances[i,j] = distance(genesP[,i+8],genesP[,j+8], genesP[,31])   
  }
}

#plot_ly(
#  x = seq(1:aux), y = seq(1:aux),
#  z = distances, type = "heatmap"
#)

rm(distances, aux, i, j, distance)

##### Integrated data #####

#Load string database
string_db = STRINGdb$new(version="10", species=9606, score_threshold=0.4, input_directory="")

head(genesP)
head(genesT)

allGenes = c(genesP$gene, genesT$gene)
allGenes = unique(allGenes)
allGenes = as.data.frame(allGenes)

allGenes = string_db$map(allGenes, "allGenes", removeUnmappedRows = TRUE)

links = string_db$get_interactions(allGenes[,2])

links = links[links$experiments!=0,c(1,2)]
head(links)

# Change the String ID for gene names
for(i in 1:nrow(links)){
  links[i,1] = allGenes[which(allGenes$STRING_id==links[i,1]),1]
  links[i,2] = allGenes[which(allGenes$STRING_id==links[i,2]),1]
}

g = graph_from_data_frame(links, directed=F)

plot(g, main = "Integrated data", vertex.color = "green", vertex.label.dist=0.3, vertex.label.degree=pi/2, vertex.size=5, vertex.label.color = "gray")

graphObj = graph_from_data_frame(links, directed=F)
layoutGraph = layout_with_fr(graphObj)

vertexCexLabel = 0.65

# Variables containing the DE genes and their logFC. "vermelho" (red in portuguese!) contains the transcipt information and "azul" , the proteomics.

vermelho = genesT[,c(1, 3)]
head(vermelho)

azul = genesP[,c(1, 3)]
head(azul)


nomes = V(graphObj)$name
nomes = as.data.frame(nomes)

nomes = merge(nomes, vermelho,  by.x = "nomes", by.y = "gene", all.x = T)
nomes = merge(nomes, azul,  by.x = "nomes", by.y = "gene", all.x = T)

head(nomes)

png(paste(folder, "histLogFCTrans.png"))
hist(nomes$logFC.x, main = "LogFC of Transcrips Distribution", col ="red", xlab = "Log FC",nc=10)
dev.off()

png(paste(folder, "histLogFCProt.png"))
hist(nomes$logFC.y, main = "LogFC of Proteomics Distribution", col ="blue", xlab = "Log FC",nc=10)
dev.off()


nomes$logFC.x = nomes$logFC.x/(max(nomes$logFC.x, na.rm=T))
nomes$logFC.x = (nomes$logFC.x-1)*(-1)

nomes$logFC.y = nomes$logFC.y/(max(nomes$logFC.y, na.rm=T))
nomes$logFC.y = (nomes$logFC.y-1)*(-1)

head(nomes)

nomes[is.na(nomes[,2]),2] = 0
nomes[is.na(nomes[,3]),3] = 0

head(nomes)

cores = c()
cores2 = c() 

for(i in 1:nrow(nomes)){
  a = nomes[i,]
  
  if(a[2]==0){
    cores = c(cores, rgb(a[3],a[3],1))
    cores2 = c(cores2, "blue")
  }else if(a[3]==0){
    cores = c(cores, rgb(1,a[2],a[2]))
    cores2 = c(cores2, "red")
  }else{
    cores = c(cores, rgb((a[2]-1)*(-1),1,(a[2]-1)*(-1)))
    cores2 = c(cores2, "green")
  }
}

head(nomes)

nomes = cbind(nomes, cores2)
nomes = cbind(nomes, cores)

nomes[1:20,]

order = c()

for(i in 1:length(V(graphObj)$name)){
  order = c(order, which(as.character(nomes[,1])==V(graphObj)$name[i]))  
}

order

head(nomes)
nomes = nomes[order,]
head(nomes)
nomes$cores = as.character(nomes$cores)

head(nomes)
nomes
situation = c()

for(i in 1:nrow(links)){
  a = as.matrix(links[i,1:2])
  aux = nomes[which(as.character(nomes[,1])==a[1]),]
  aux2 = nomes[which(as.character(nomes[,1])==a[2]),]
  
  if(aux[2]==0 && aux2[2]==0){
    situation = c(situation,"blue")
  } else if(aux[3]==0 && aux2[3]==0){
    situation = c(situation,"red")
  } else {
    situation = c(situation,"green")
  }
}

table(situation)
situation

# Plot with colors and all conections

#png(paste(folder, "completeNetwork.png", sep=""))
plot(graphObj, main = "Integrated transcript and proteomic network", vertex.color = nomes$cores, vertex.size=5, 
     edge.color = situation, vertex.label.dist=0.25, vertex.label.degree=pi/2,
     vertex.label.cex=vertexCexLabel, vertex.label.color = "gray", layout = layoutGraph)
legend(0.65,-1, c("Proteomics","Transcripts",  "Both"), fill = c("blue", "red", "green"), cex = 0.75, title = "Edges")
if(zeroT>0.8){
  col.labels<-c("Min", "", "","","Max/Zero")  
}else if(zeroT>0.6){
  col.labels<-c("Min", "", "","Zero","Max")
}else if(zeroT>0.4){
  col.labels<-c("Min", "", "Zero","","Max")
}else if(zeroT>0.2){
  col.labels<-c("Min", "Zero", "","","Max")
}else{
  col.labels<-c("Min/Zero", "", "","","Max")
}
testcol<-color.gradient(1,c(1,0),c(1,0),nslices=100)
# color.legend(-1,-0.85,-0.85,-1.25,col.labels,testcol,gradient="y")
# text(-1,-0.75,"Transcripts", cex=0.6)
# text(-1,-0.8,"Fold Change", cex=0.6)
color.legend(1,1,1.15,0.6,col.labels,testcol,gradient="y", cex = 0.6)
text(1,1.1,"Transcripts", cex=0.6)
text(1,1.05,"Fold Change", cex=0.6)

if(zeroP>0.8){
  col.labels<-c("Min", "", "","","Max/Zero")  
}else if(zeroP>0.6){
  col.labels<-c("Min", "", "","Zero","Max")
}else if(zeroP>0.4){
  col.labels<-c("Min", "", "Zero","","Max")
}else if(zeroP>0.2){
  col.labels<-c("Min", "Zero", "","","Max")
}else{
  col.labels<-c("Min/Zero", "", "","","Max")
}
testcol<-color.gradient(c(1,0),c(1,0),1,nslices=100)
color.legend(-1,1.00,-0.85,0.6,col.labels,testcol,gradient="y", cex=0.6)
text(-1,1.10,"Proteomics", cex=0.6)
text(-1,1.05,"Fold Change", cex=0.6)
#dev.off()

##### Second plot (same plot, but with only the transcripts connections) #####

sit = situation
sit[sit=="blue"] = "white"
sit[sit=="green"] = "white"

nom = nomes
dummyIdx = which(nom$cores2=="blue")
nom$cores2[dummyIdx] = "white"
nom$cores[dummyIdx] = "#FFFFFF"

vertexFrameColor = rep("black", length = dim(nom)[1], mode = "character")
vertexFrameColor[dummyIdx] = "white"

labelColor = rep("gray", length = dim(nom)[1], mode = "character")
labelColor[dummyIdx] = "white"

plot(graphObj, main = "Transcript network", vertex.color = nom$cores, vertex.frame.color = vertexFrameColor, 
     vertex.size=5, edge.color = sit, vertex.label.dist=0.25, vertex.label.cex=vertexCexLabel,
     vertex.label.degree=pi/2, vertex.label.color = labelColor, layout = layoutGraph)
# legend("bottomright", c("Proteomics","Transcripts",  "Both"), fill = c("blue", "red", "green"))

# oether way deleteing esge
graphObjT = graphObj
delete_edges(graphObjT, which((situation=="blue")&(situation=="green")))
delete_vertices(graphObjT, nomes[nomes$cores2=="red",1])
# sit = situation[situation=="red"]
layoutGraphT = layoutGraph[dummyIdx,]
plot(graphObjT, main = "Transcript network", vertex.color = nom$cores, vertex.frame.color = vertexFrameColor, 
     vertex.size=5, edge.color = "red", vertex.label.dist=0.25, vertex.label.cex=vertexCexLabel,
     vertex.label.degree=pi/2, vertex.label.color = "grey", layout = layoutGraph)

##### Third plot (same plot, but with only the proteomics connections) #####

sit = situation

sit[sit=="red"] = "white"
sit[sit=="green"] = "white"

nom = nomes
dummyIdx = which(nom$cores2=="red")
nom$cores2[dummyIdx] = "white"
nom$cores[dummyIdx] = "#FFFFFF"

vertexFrameColor = rep("black", length = dim(nom)[1], mode = "character")
vertexFrameColor[dummyIdx] = "white"

labelColor = rep("gray", length = dim(nom)[1], mode = "character")
labelColor[dummyIdx] = "white"

plot(graphObj, main = "Proteomic network", vertex.color = nom$cores, vertex.frame.color = vertexFrameColor,
     vertex.size=5, edge.color = sit, vertex.label.dist=0.25, vertex.label.degree=pi/2,
     vertex.label.cex=vertexCexLabel, vertex.label.color = labelColor, layout = layoutGraph)
# legend("bottomright", c("Proteomics","Transcripts",  "Both"), fill = c("blue", "red", "green"))

##### GO information and fourth plot  - Under construction #####

head(genesT)

annotations = string_db$get_annotations()

go =  string_db$get_enrichment(genesT[,27])

head(annotations)
head(go)

genesT = merge(genesT, annotations, by.x = "STRING_id", by.y = "STRING_id")

head(genesT)

genesT = merge(genesT, go, by.x = "term_id", by.y = "term_id")

head(genesT)

unique(genesT$term_description)

# Choose biological function here
# bio = "immune response"
# bio = "response to reactive oxygen species"
# bio = "response to interleukin-1"
# bio = "leukocyte migration"
bio = "extracellular matrix organization"
# bio = "inflammatory response"
#bio = "response to cytokine"

head(nomes)

aux = genesT[genesT$term_description==bio,3]

aux = as.character(aux)

for(i in 1:length(aux)){
  aux[i] = which(V(graphObj)$name==aux[i])
}

aux = as.numeric(aux)

aux2 = rep("gray", length(V(graphObj)$name))

aux2[aux] = "purple"

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(graphObj, main = "Gene Network", vertex.color = nomes$cores, vertex.size=5, edge.color = situation, 
     vertex.label.dist=0.25, vertex.label.degree=pi/2, vertex.label.color=aux2, layout = layoutGraph,
     vertex.label.cex=vertexCexLabel)
legend("bottomright", inset=c(-0.2,0), c("Proteomics","Transcripts",  "Both", bio), fill = c("blue", "red", "green", "purple"), cex = 0.75)



