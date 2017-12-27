head(dados)

pedig = with(dados, pedigree(id=dados$IID, dadid=dados$PAT, momid=dados$MAT, 
                             sex=(dados$SEX), famid=dados$FID, missid=0))
kmat = kinship(pedig)

a = pedig['14']
plot(a)


exemplo = cbind(c(1:9), c(0,0,0,2,2,2,0,3,7), c(0,0,0,1,1,1,0,4,6),c(1,0,0,1,0,1,0,1,0), rep(1,9))
exemplo = as.data.frame(exemplo)
colnames(exemplo) = c("ID","PAT","MAT", "SEX", "FID")
exemplo


pedig = pedigree(id=exemplo$ID, dadid=exemplo$PAT, momid=exemplo$MAT, 
                               sex=exemplo$SEX, famid=exemplo$FID, missid=0)

par(bg=NA)
plot(pedig['1'], main = "Family Pedigree")

kmat=kinship(pedig)
kmat
