source("kclust.R")
set.seed(0)

x = rbind(matrix(rnorm(2*100,sd=0.2),ncol=2),
scale(matrix(rnorm(2*100,sd=0.3),ncol=2),cent=-c(1,1),scal=F),
scale(matrix(rnorm(2*100,sd=0.2),ncol=2),cent=-c(0,1),scal=F))

k = 3
cent.init = rbind(c(0.5,1),c(1,0),c(0,0.5))

cols = c("red","darkgreen","blue")
plot(x)
points(cent.init,pch=19,cex=2,col=cols)

km1 = kclust(x,centers=cent.init,alg="kmeans")
km2 = kmeans(x,centers=cent.init,alg="Lloyd")
sum(km1$cluster!=km2$cluster)

plot(x,col=cols[km1$cluster])
points(km1$centers,pch=19,cex=2,col=cols)

cent.old = cent.init

plot(x)
points(cent.old,pch=19,cex=2,col=cols)

par(ask=TRUE)

for (i in 1:km1$iter) {
  # Plot the new clusters
  plot(x,col=cols[km1$cluster.history[i,]],main=paste("Iteration",i))
  points(cent.old,pch=19,cex=2,col=cols) 

  # Plot the new centers
  plot(x,col=cols[km1$cluster.history[i,]],main=paste("Iteration",i))
  cent.new = centfromclust(x,km1$cluster.history[i,],k,alg="kmeans")
  points(cent.new,pch=19,cex=2,col=cols)

  cent.old = cent.new
}




