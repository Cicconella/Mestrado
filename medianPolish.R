
a = matrix(round(runif(50, 1, 15)), ncol=10)
a[1,5] = 20

sub = a
sub

mediana1 = mediana2 = 5

while(sum(abs(round(mediana1, digits = 2)))>0 & sum(abs(round(mediana2, digits = 2)))>0){
  mediana1 = apply(sub,1,median)
  mediana1
  sub = sub - mediana1
  sub
  
  mediana2 = apply(sub,2,median)
  mediana2
  
  for(i in 1:10){
    sub[,i] = sub[,i]-mediana2[i]
  }
  sub  
}

final = a - sub
