
a = matrix(round(runif(50, 1, 15)), ncol=10)
a[1,5] = 20

a = matrix(c(1,2,3,3,4,5,5,6,7), nrow=3, byrow = T)
medpolish(a)

matrix(c(1,2,3), nrow = 3, ncol = 3)

medpolish

sub = a
sub

mediana1 = mediana2 = 0

while(sum(abs(round(mediana1, digits = 2)))>0 & sum(abs(round(mediana2, digits = 2)))>0){
  mediana1 = apply(sub,1,median)
  mediana1
  sub = sub - mediana1
  sub
  
  mediana2 = c(apply(sub,2,median), median(mediana1))
  mediana2
  
  for(i in 1:10){
    sub[,i] = sub[,i]-mediana2[i]
  }
  sub  
}

final = a - sub
x=a

z <- as.matrix(x)
nr <- nrow(z)
nc <- ncol(z)
t <- 0
r <- numeric(nr)
c <- numeric(nc)
oldsum <- 0
maxiter = 10L
for (iter in 1L:maxiter) {
  rdelta <- apply(z, 1L, median, na.rm = F)
  z <- z - matrix(rdelta, nrow = nr, ncol = nc)
  r <- r + rdelta
  delta <- median(c, na.rm = F)
  c <- c - delta
  t <- t + delta
  cdelta <- apply(z, 2L, median, na.rm = T)
  z <- z - matrix(cdelta, nrow = nr, ncol = nc, byrow = TRUE)
  c <- c + cdelta
  delta <- median(r, na.rm = T)
  r <- r - delta
  t <- t + delta
  newsum <- sum(abs(z), na.rm = T)
  converged <- newsum == 0 || abs(newsum - oldsum) < eps * 
    newsum
  if (converged) 
    break
  oldsum <- newsum
  if (trace.iter) 
    cat(iter, ": ", newsum, "\n", sep = "")
}
if (converged) {
  if (trace.iter) 
    cat("Final: ", newsum, "\n", sep = "")
}
else warning(sprintf(ngettext(maxiter, "medpolish() did not converge in %d iteration", 
                              "medpolish() did not converge in %d iterations"), maxiter), 
             domain = NA)
names(r) <- rownames(z)
names(c) <- colnames(z)
ans <- list(overall = t, row = r, col = c, residuals = z, 
            name = deparse(substitute(x)))
class(ans) <- "medpolish"
ans






