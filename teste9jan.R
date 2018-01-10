for(i in 1:nrow(cnv)){
  
  cn = cnv[i,-c(1:3)]
  
  cn = cbind(colnames(cn),as.numeric(as.character(cn)))
  cn = as.data.frame(cn)
  head(cn)
  summary(cn)
  
  class(cn[,2])
  
  dadoscn = merge(dados, cn, by.x="cel", by.y="V1", all.x = T)
  
  head(dadoscn)
  
  class(cn$V2)
  table(cn$V2)
  
  aux = table(cn$V2)
  print(i)    
  print(aux)

  if(length(which(aux<5))!=0){
    print("tem que tirar")
  }

}
 