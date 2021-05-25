
library(car)
library(pls)
library(parcor)
library(lars)

setwd("Adresse de votre fichier")
cookies = read.table("cookies.txt", sep=" ",header=TRUE)

n = nrow(cookies)
p=ncol(cookies)


Y=as.matrix(cookies[1])
X=as.matrix(cookies[-1])

plot(X[1,1],Y[1],xlab="Spectre",ylab="Teneur en sucre",xlim=c(0,1),ylim=c(8,23))
for (j in 1:10){
for (i in 1:200) points(X[4*j,i],Y[j],col=j)
}
title("Observation de 10 cookies")
Cor=rep(1,10)
for(i in 1:10) {Cor[i]=cor(X[,i],X[,i+1])}



# Estimation PLS.
RegPLS = plsr(Y ~ X, ncomp=30, scale=TRUE, validation='LOO')
plot(MSEP(RegPLS, estimate='CV'),xlab="Nombre de composantes")
RegPLSopt = plsr(Y ~ X, ncomp=21, scale=TRUE)
title("Détermination de k")

coef(RegPLSopt)

# PRESS PLS
PredPLS = rep(0, n)
for (k in 1:n){
  RegPLSk = plsr(Y[-k] ~ X[-k,], ncomp=21, scale=TRUE)	
  
  X2 = X
  for (i in 1:700){
  XCol = X2[,i]
  X2[,i] = (XCol-mean(XCol[-k]))/sd(XCol[-k])
  }
  PredPLS[k] = coef(RegPLSk)%*%c(X2[k,]) + mean(Y[-k])
}
##### Press
PRESSPLS = sum((Y - PredPLS)^2)


##### Graph PLS
plot(1:n,Y, type="o", lwd=1.5,ylim=c(8,25),ylab=("Teneur en sucre"))
lines(1:n, PredPLS, type="o", col="red", lty=2)
legend("topleft", c("Obs.", "PLS"), lty=1, col=c("black", "red"),bty="n")
title("Estimation PLS")


# Estimation Ridge
RCVRidge = ridge.cv(X, Y, lambda=seq(0.1, 50, 0.1), k=nrow(X), plot.it=TRUE)
abline(v=7.2,col="red",lty=1)
lambda = RCVRidge$lambda.opt
RegRidge = lm.ridge(Y ~ X, lambda=lambda)

# PRESS Ridge
PredRidge = rep(0, n)
for (k in 1:n){
  RegRidgek = lm.ridge(Y[-k] ~ X[-k,], lambda=lambda)	
  PredRidge[k] = coef(RegRidgek)%*%c(1, X[k,])
}
PRESSRidge = sum((Y - PredRidge)^2)

## Graph Ridge
plot(1:n,Y, type="o", lwd=1.5,ylim=c(7,25),ylab=("Teneur en sucre"))
lines(1:n, PredRidge, type="o", col="red", lty=2)
legend("bottomleft", c("Obs.", "Ridge"), lty=1, col=c("black", "red"),bty="n")
title("Estimation Ridge")


# Estimation Lasso
RCVLasso = cv.lars(X, Y, K=nrow(X), type="lasso",use.Gram=F)
abline(v=1/3,col="red",lty=1)
delta = RCVLasso$index[which.min(RCVLasso$cv)]
RegLasso = lars(X, Y,use.Gram=F)

# PRESS Lasso
PredLasso = rep(0, n)
for (k in 1:n){
  X2 = X
  for (i in 1:700){
  XCol = X2[,i]
  X2[,i] = XCol-mean(XCol[-k])
  }
  PredLasso[k] = coef(RegLassok, s=delta, mode="fraction")%*%X2[k,] + mean(Y[-k])
}
PRESSLasso = sum((Y - PredLasso)^2)

plot(1:n, Y, type="o", lwd=1.5)
lines(1:n, PredLasso, type="o", col="red", lty=2)
legend("topleft", c("Obs.","Lasso"), lty=1, col=c("black", "red"),bty="n")
title("Estimation Lasso")


#####################


plot(1:n, Y, type="o", lwd=2)
lines(1:n, PredPLS, type="o", col="blue", lty=2)
lines(1:n, PredRidge, type="o", col="red", lty=2)
lines(1:n, PredLasso, type="o", col="chartreuse4", lty=2)
legend("topleft", c("Obs.", "PLS", "Ridge", "Lasso"), lty=1, col=c("black", "blue", "red", "chartreuse4"))
title("Résumé des estimations")


##################### Nouveau jeu de données

newcook=read.table("nouveaux.txt",sep=" ",header=TRUE)
T=as.matrix(newcook[-1])

##################### PLS

PredPLS2 = rep(0, 32)
for (k in 1:32){
  X2 = T
  for (i in 1:700){
    XCol = X2[,i]
    X2[,i] = (XCol-mean(XCol[-k]))/sd(XCol[-k])
  }
  PredPLS2[k] = coef(RegPLS)%*%X2[k,] + mean(Y[-k])
}


################# LASSO

PredLasso2=rep(0,32)
for (k in 1:32) {
  X2 = T
  for (i in 1:700){
  XCol = X2[,i]
  X2[,i] = XCol-mean(XCol[-k])
  }
  PredLasso2[k] = coef(RegLasso, s=delta, mode="fraction")%*%X2[k,] + mean(Y[-k])
}


################# RIDGE

PredRidge2=rep(0,32)
for (k in 1:32){
  PredRidge2[k] = coef(RegRidge)%*%c(1,T[k,])
}


plot(1:32, newcook[,1], type="o", lwd=1.5,xlab="n°",ylab="Teneur en sucre",ylim=c(8,25))
lines(1:32, PredLasso2, type="o", col="red", lty=2)
legend("bottomleft", c("Obs.","Lasso"), lty=1, col=c("black","red"),bty="n")
title("Prévision Lasso")

#################

sum((newcook[,1]-PredPLS2)^2)
sum((newcook[,1]-PredRidge2)^2)
sum((newcook[,1]-PredLasso2)^2)
