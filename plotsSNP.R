
a = read.table("intensity-snp")

head(a)
dim(a)

a[1:2,1:10]

plot(t(a[1:2,-1]), xlab = "Allele A", ylab = "Allele B", pch = 16, col = rgb(0,0,1, alpha = 0.5), 
     main = "Allele intensities of\nSNP_A-4265735", cex.lab = 1.5, cex.main=1.5)

r = a[1,-1]+a[2,-1]
r[1:5]
r = as.numeric(r)

t = (atan(a[2,-1]/a[1,-1]))/(pi/2)
t[1:5]
t = as.numeric(t)

plot(r,t, xlab = "R", ylab = expression(theta), pch = 16, col = rgb(0,0,1, alpha = 0.5), 
     main = ("Polar coordinates of\nSNP_A-4265735 intensities"), cex.lab = 1.5, cex.main = 1.5)

a = read.table("/home/cicconella/DadosMestrado/lrr-baf-snp")

dim(a)
a[1:10]

lrr = a[seq(4,2243, by=2)]
baf = a[seq(5,2243, by=2)]

plot(as.numeric(lrr), as.numeric(baf), pch = 16,  col = rgb(0,0,1, alpha = 0.5), xlab = "Log R Ratio", ylab = "BAF",
     main = "HMM parameters for SNP_A-4265735")
