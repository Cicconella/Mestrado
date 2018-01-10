

samples = rbind(c(123,124,125),
                c(124,0,0),
                c(125,0,0),
                c(122,124,125))
colnames(samples) = c("ID","FA", "MO")
samples = as.data.frame(samples)
samples

cnv = cbind(c(1,1,2,2),
            c(2,2,2,2),
            c(3,2,NA,2))
colnames(cnv) = c("CN1","CN2","CN3")
cnv = as.data.frame(cnv)
cnv

c=1
for(c in 1:3){
  cn = as.numeric(cnv[,c])
  cn = cbind(samples,cn)
  
  comb = expand.grid(c(0:4),c(0:4),c(0:4))
  comb = data.frame(cbind(comb, rep(0,125)))
  comb[,4] = as.numeric(as.character(comb[,4]))
    
  colnames(comb) = c("P1", "P2", "OF", "CN")
    
  head(comb)
  head(cn)
  i=4
  for(i in 1:nrow(cn)){
      cn[i,]
      
      if(is.na(as.character(cn[i,"cn"]))){
        print("no cnv information")
        next
      }
      
      if(cn[i,"FA"]!=0 & cn[i,"MO"]!=0){
        if(is.na(as.character(cn[which(cn[,1]==cn[i,"FA"]),"cn"]))
           | is.na(as.character(cn[which(cn[,1]==cn[i,"MO"]),"cn"]))){
          next
        }else{
          a = t(as.matrix(c(as.numeric(as.character(cn[which(cn[,1]==cn[i,"FA"]),"cn"])),
                            as.numeric(as.character(cn[which(cn[,1]==cn[i,"MO"]),"cn"])),
                            as.numeric(as.character(cn[i,"cn"])))))
          colnames(a) = colnames(comb)[-4]
          comb[which(apply(comb, 1, function(x) identical(x[1:3], a[1,]))),4] = comb[which(apply(comb, 1, function(x) identical(x[1:3], a[1,]))),4]+1
        }
      }else{
        print("no parents")
        next
      }
  }
    comb
    
    for(i in 1:75){
      a = t(as.matrix(as.numeric(c(rev(comb[i, 1:2]),comb[i,3]))))
      colnames(a) = c("P1", "P2","OF")
      bla = which(apply(comb, 1, function(x) identical(x[1:3], a[1,])))
      print(bla)
      print(i)
      if(bla!=i){
        comb[i,4] = comb[i,4]+comb[bla,4]
        comb = comb[-bla,]
      }
    }
    
    dim(comb)
    
    
    if(c == 1){
      trios = comb
    }else{
      trios = cbind(trios,comb[,4])  
    }
  
}

samples
cnv
trios
  