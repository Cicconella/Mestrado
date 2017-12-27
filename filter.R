cont_0 = function(x){
  length(which(x[-c(1,2,3)]==0))
}

cont_1 = function(x){
  length(which(x[-c(1,2,3)]==1))
}

cont_2 = function(x){
  length(which(x[-c(1,2,3)]==2))
}

cont_3 = function(x){
  length(which(x[-c(1,2,3)]==3))
}
cont_4 = function(x){
  length(which(x[-c(1,2,3)]==4))
}


grupos = function(a){
  t = sum(a)
  g = 5
  
  if(a[1]<t*mcf)
    g = g-1
  if(a[2]<t*mcf)
    g = g-1
  if(a[3]<t*mcf)
    g = g-1
  if(a[4]<t*mcf)
    g = g-1
  if(a[5]<t*mcf)
    g = g-1
  
  return(g)
}


cnv = cnvA
cnv




mcf = 0.02


mut_0 = apply(cnv,1,cont_0)
mut_1 = apply(cnv,1,cont_1)
mut_2 = apply(cnv,1,cont_2)
mut_3 = apply(cnv,1,cont_3)
mut_4 = apply(cnv,1,cont_4)

mutations = cbind(mut_0,mut_1,mut_2,mut_3,mut_4)

dim(cnv)
dim(mutations)
  
groupCNV = apply(mutations, 1, grupos)
  
summary(as.factor(groupCNV))

exclude = which(groupCNV==1)
cnv = cnv[-exclude,]
dim(cnv)
cnv
  
