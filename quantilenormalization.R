a = round(runif(10, min = 0, max=10))
#[1]  0  4  3  8  5  1  2  6  5 10

b = round(runif(10, min = 5, max=10))
#[1] 10  8  9  9  8  9  9  9  7  9

plot(a, b, pch="x", xlim=c(0,10), ylim=c(0,10), xlab = "Array A", 
     ylab = "Array B", main = "Raw array values")
lines(0:10,0:10)


qqplot(c,d, pch="x", xlim=c(0,10), ylim=c(0,10), xlab = "Sorted array A", 
     ylab = "Sorted array B", main = "Quantile-quantile plot")
lines(0:10,0:10)

e = cbind(sort(a),sort(b))
f = apply(e, 1, mean)

qqplot(f[order(a)],f[order(b)], pch="x", xlim=c(0,10), ylim=c(0,10),
       xlab = "Normalized array A", ylab = "Normalized array B",
       main = "Quantile-quantile plot")
lines(0:10,0:10)

f[order(a)]
# [1]  3.5  7.0  7.0  5.0  4.5  6.5  8.5  7.5  6.0 10.0
f[order(b)]
# [1]  8.5  4.5  6.5  5.0  6.0  7.0  7.0  7.5 10.0  3.5