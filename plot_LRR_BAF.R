plot(1:10, (-4:5)^2, main = "Parabola Points", xlab = "xlab")
mtext("10 of them")
for(s in 1:4)
  mtext(paste("mtext(..., line= -1, {side, col, font} = ", s,
              ", cex = ", (1+s)/2, ")"), line = -1,
        side = s, col = s, font = s, cex = (1+s)/2)
mtext("mtext(..., line= -2)", line = -2)
mtext("mtext(..., line= -2, adj = 0)", line = -2, adj = 0)
##--- log axis :
plot(1:10, exp(1:10), log = "y", main = "log =\"y\"", xlab = "xlab")
for(s in 1:4) mtext(paste("mtext(...,side=", s ,")"), side = s)


x = c(sample(c(runif(10,0.45,0.55),rep(0,10),rep(1,10))),
      sample(c(runif(10,0.45,0.55),runif(10,0.22,0.28), runif(10,0.72,0.78),rep(0,10),rep(1,10))),
      sample(c(runif(10,0.45,0.55),rep(0,10),rep(1,10))),
      sample(c(runif(15,0.3,0.36),runif(15,0.63,0.69),rep(0,10),rep(1,10))),
      sample(c(runif(10,0.45,0.55),rep(0,10),rep(1,10))),
      sample(c(rep(0,25),rep(1,25))),
      sample(c(runif(10,0.45,0.55),rep(0,10),rep(1,10))),
      sample(c(runif(40,0,1),rep(0,5),rep(1,5))),
      sample(c(runif(10,0.45,0.55),rep(0,10),rep(1,10))),
      sample(c(rep(0,25),rep(1,25))),
      sample(c(runif(10,0.45,0.55),rep(0,10),rep(1,10))),
      sample(c(runif(15,0.3,0.36),runif(15,0.63,0.69),rep(0,10),rep(1,10))),
      sample(c(runif(10,0.45,0.55),rep(0,10),rep(1,10))),
      sample(c(runif(15,0.38,0.4),runif(15,0.55,0.61),rep(0,10),rep(1,10))),
      sample(c(runif(10,0.45,0.55),rep(0,10),rep(1,10))))


plot(NULL, yaxt="n", xaxt="n", xlim=c(0,length(x)), ylim=c(0,1), 
     xlab="Position", ylab="BAF", cex.lab=1.4)
axis(2,c(0,0.5,1))
abline(v=seq(30.5,80.5),col="lightblue")
abline(v=seq(110.5,160.5),col="lightblue")
abline(v=seq(190.5,240.5),col="lightblue")
abline(v=seq(270.5,320.5),col="lightblue")
abline(v=seq(350.5,400.5),col="lightblue")
abline(v=seq(430.5,480.5),col="lightblue")
abline(v=seq(510.5,560.5),col="lightblue")
lines(x=seq(1,length(x)),y=rep(0.5,length(x)))
points(x,pch=20)



x = c(runif(30,-0.05,0.05),
      runif(50,0.2,0.3),
      runif(30,-0.05,0.05),
      runif(50,0.1,0.2),
      runif(30,-0.05,0.05),
      runif(50,0.1,0.2)*(-1),
      runif(30,-0.05,0.05),
      runif(50,0.2,0.3)*(-1),
      runif(30,-0.05,0.05),
      runif(50,-0.05,0.05),
      runif(30,-0.05,0.05),
      runif(50,0,0.1)*(-1),
      runif(30,-0.05,0.05),
      runif(50,-0.05,0.1),
      runif(30,-0.05,0.05))


plot(NULL, yaxt="n", xaxt="n", xlim=c(0,length(x)), ylim=c(-0.4,0.4), 
     cex.lab=1.4, ylab="LRR", xlab="")
axis(2,c(0))
abline(v=seq(30.5,80.5),col="lightblue")
abline(v=seq(110.5,160.5),col="lightblue")
abline(v=seq(190.5,240.5),col="lightblue")
abline(v=seq(270.5,320.5),col="lightblue")
abline(v=seq(350.5,400.5),col="lightblue")
abline(v=seq(430.5,480.5),col="lightblue")
abline(v=seq(510.5,560.5),col="lightblue")
lines(x=seq(1,length(x)),y=rep(0,length(x)))
points(x,pch=20)

