a = round(runif(5, min = 0, max=10))
b = round(runif(5, min = 5, max=10))

plot(a,b, pch="x", xlim=c(0,10), ylim=c(0,10))
lines(0:10,0:10)

c = sort(a)
d = sort(b)

plot(c,d, pch="x", xlim=c(0,10), ylim=c(0,10))
lines(0:10,0:10)

qqplot(a,b, pch="x", xlim=c(0,10), ylim=c(0,10))
lines(0:10,0:10)

e = cbind(c,d)

medias = apply(e,1,mean)

f = medias[order(a)]
g = medias[order(b)]

plot(f, g, pch="x")
qqplot(f,g)

a
c
f[,1]

b
d
f[,2]

ordem1 = c(5,9,10,9,3,1,7,4,6,2)
ordem2 = c(2,6,5,3,7,9,4,8,1,10)

g = f[ordem1,1] 
h = f[ordem2,2] 

plot(g,h)
qqplot(g,h)
