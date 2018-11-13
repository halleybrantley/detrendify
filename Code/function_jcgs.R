# This file includes functions for nonparametric quantile estimation (JCGS, 2010)

hry1 <- function(k,n,yux){
out <- vector(length=n-k+1,mode="numeric")
for(i in 1:(n-k+1)){
out[i] <- median(yux[i:(i+k-1)])}
out
}


hry2 <- function(k,p,n,yuy){
out <- vector(length=n-k+1,mode="numeric")
for(i in 1:(n-k+1)){
out[i] <- quantile(yuy[i:(i+k-1)],prob=p)}
out
}

aresult <- function(xrand, n, h, kemx, kemy){
       lenx <- length(xrand)
       out <- vector(length = lenx, mode = "numeric")
       for(i in 1:lenx) {
               out[i] <- bigfun1(n, xrand[i], h, kemx, kemy)
       }
       out
}


smoo2 <- function(n, h, x, i, kemx){
       vec1 <- rep(1, n)
       xvec1 <- x * vec1
       xvec2 <- kemx - xvec1
       out <- as.vector(xvec2^i %*% dnorm(xvec2/h)/h)
       out
}

bigfun1 <- function(n, x, h, kemx, kemy){
       out <- (smoo2(n, h, x, 2, kemx) * tfun1(n, 0, x, h, kemx, kemy) - smoo2(
               n, h, x, 1, kemx) * tfun1(n, 1, x, h, kemx, kemy))/(smoo2(n, h,
               x, 2, kemx) * smoo2(n, h, x, 0, kemx) - smoo2(n, h, x, 1, kemx)^
               2 + n^-5)
       out
}

tfun1 <- function(n, i, x, h, kemx, kemy){
       vec1 <- rep(1, n)
       xvec1 <- x * vec1
       xvec2 <- kemx - xvec1
       yvec1 <- kemy * xvec2^i
       out <- as.vector(yvec1 %*% dnorm(xvec2/h)/h)
       out
}


# an approximate rho function
rho.m <- function(x,tau,cut){
   tmp1 <- tau*x[x>=cut]-0.5*cut*tau
   tmp2 <- 0.5*tau*x[0<=x & x<cut]^2/cut
   tmp3 <- 0.5*(1-tau)*x[(-cut)<=x & x<0]^2/cut
   tmp4 <- (1-tau)*abs(x[x<(-cut)])-0.5*cut*(1-tau)
   rho.f <- c(tmp4,tmp3,tmp2,tmp1)
   rho.f
}
# the derivative of rho.m
dm.rho <- function(x,tau,cut){
   x <- pmax(x,-cut)
   x <- pmin(x,cut)
   x <- x
   x[x>0] <- tau*x[x>0]
   x[x<=0] <- (1-tau)*x[x<=0]
   x
}
# NOTE: derivative of rho.m() should be dm.rho()/cut
#       as another scaled version, we can consider
#       dm.rho()/(tau/cut)



# k-NN method for initial fit with boundary correction
# symmetric case
int.fit.sy <- function(k,n,p,ry){
   qp <- NULL
   v.st <- median(ry[1:floor(k/10)])
   v.end <- median(ry[n:(n-floor(k/10))])
   #ny <- c(v.st-ry[k:1],ry,v.end+ry[(n-(k-1)):n])
   ny <- c(v.st-(ry[k:1]-v.st),ry,v.end+(ry[(n-(k-1)):n]-v.end))
   for(i in c(1:n)){
       i.y <- quantile(ny[i:(i+2*k-1)],prob=p)
       qp <-c(qp,i.y)
   }
qp
}



int.fit.fig <- function(k,n,p,ry){
   qp <- NULL
   v.st <- median(ry[1:floor(k/10)])
   v.end <- median(ry[n:(n-floor(k/10))])
   ny <- c(v.st-(ry[(k+1):2]-v.st),ry,v.end+(v.end-ry[(n-1):(n-k)]))
   for(i in c(1:n)){
       i.y <- quantile(ny[i:(i+2*k-1)],prob=p)
       qp <-c(qp,i.y)
   }
qp
}


# get an initial fit in 2-d case
d2int.fit.sy <- function(x,y,z,k,p){
   n <- length(x)
   m <- length(y)
   z.exp <- matrix(0,(n+2*k),(m+2*k))
   nn <- nrow(z.exp)
   mm <- ncol(z.exp)
   z.exp[(k+1):(nn-k),(k+1):(mm-k)] <- z
   tmp <- rep(1,k)%*%z[1:k,]
   tmp <- rep(tmp/k,k)
   tmp.z <- matrix(tmp,nrow=k,byrow=T)
   z.exp[(1:k),(k+1):(mm-k)] <- tmp.z-(z[k:1,]-tmp.z)
   tmp <- z[,1:k]%*%rep(1,k)
   tmp <- rep(tmp/k,k)
   tmp.z <- matrix(tmp,ncol=k)
   z.exp[(k+1):(nn-k),(1:k)] <- tmp.z-(z[,(k:1)]-tmp.z)
   al <- nn-k+1
   bl <- n-k+1
   tmp <- rep(1,k)%*%z[bl:n,]
   tmp <- rep(tmp/k,k)
   tmp.z <- matrix(tmp,nrow=k,byrow=T)
   z.exp[al:nn,(k+1):(mm-k)] <- tmp.z+(z[(bl:n),]-tmp.z)
   cl <- mm-k+1
   dl <- m-k+1
   tmp <- z[,bl:m]%*%rep(1,k)
   tmp <- rep(tmp/k,k)
   tmp.z <- matrix(tmp,ncol=k)
   z.exp[(k+1):(nn-k),cl:mm] <- tmp.z+(z[,(dl:m)]-tmp.z)
   z.exp[1:k,1:k] <- 0.5*(z.exp[1:k,(2*k):(k+1)]+z.exp[(2*k):(k+1),1:k])
   z.exp[1:k,cl:mm] <- 0.5*(z.exp[1:k,(mm-k):(cl-k)]+z.exp[(2*k):(k+1),cl:mm])
   z.exp[al:nn,1:k] <- 0.5*(z.exp[(nn-k):(al-k),1:k]+z.exp[al:nn,(2*k):(k+1)])
   z.exp[al:nn,cl:mm] <- 0.5*(z.exp[(nn-k):(al-k),(cl:mm)]+z.exp[al:nn,(mm-k):(cl-k)])
   zint <- matrix(0,n,m)
   for(i in 1:n){
       for(j in 1:m){
           tmp <- as.vector(z.exp[i:(i+2*k-1),j:(j+2*k-1)])
           zint[i,j] <- quantile(tmp,prob=p)
       }
   }
   return(zint)
}



# periodic case
int.fit.pe <- function(k,n,p,ry){
   qp <- NULL
   ny <- c(ry[1]-ry[k:1],ry,ry[n]+ry[(n-(k-1)):n])
   for(i in c(1:n)){
       i.y <- quantile(ny[i:(i+2*k-1)],prob=p)
       qp <-c(qp,i.y)
   }
qp
}


# quantile smoothing spline
qreq1d.sreg<-function(x,y,finit,p=0.5,cutoff=0.01) {
 n<-length(y)
 fnew <- finit
 pynew <- y
 err<-1
 niter<-0
 while ((err>0.01)&(niter<30)) {
   fold <- fnew
   pyold <- pynew
   resid <- pyold-fold
   pynew <- fold + dm.rho(resid,p,cutoff)
   temp <- sreg(x,pynew)
   fnew <- temp$predicted$y
   err<-mean(abs(fnew-fold)/pmax(abs(fnew),0.01))
   niter<-niter+1
#    cat(niter, err, "\n")
#    plot(x,pynew,pch=".")
#    lines(x,fold)
 }
 return(fnew)
}

# it is used for simulation of 1-d case in the paper.
qreq1d.sreg.fig <- function(x,y,finit,tau,cut) {
       n<-length(y)
       fnew <- finit
       pynew <- y
       err<-1
       niter<-0
       while ((err>0.01)&(niter<30)) {
       fold <- fnew
       pyold <- pynew
       resid <- pyold-fold
       pynew <- fold + dm.rho(resid,tau,cut*10*x+0.1)/0.2
       temp <- sreg(x,pynew)
       fnew <- temp$fitted.values
       err <- mean((fnew-fold)^2)
       niter<-niter+1
#    cat(niter, err, "\n")
#    plot(x,pynew,pch=".")
#    lines(x,fold)
       }
       return(fnew)
}



# quantile wavelet shrinkage
qreq1d.wavelet<-function(y,finit,p=0.5,cutoff=0.01) {
 # assume y is equi-spaced
 n<-length(y)
 fnew <- finit
 pynew <- y
 err<-1
 niter<-0
 while ((err>0.01)&(niter<30)) {
   fold <- fnew
   pyold <- pynew
   resid <- pyold-fold
   pynew <- fold + dm.rho(resid,p,cutoff)
   fnew <- wr(threshold(wd(pynew,filter.number=8,family=c("DaubLeAsymm")),policy="universal"))
   err<-mean(abs(fnew-fold)/pmax(abs(fnew),0.01))
   niter<-niter+1
   cat(niter, err, "\n")
#    plot(x,pynew,pch=".")
#    lines(x,fold)
 }
 return(fnew)
}


med2dknn<-function(x,y,z,k) {
 n<-length(x)
 medz<-z
 for (i in 1:n) {
   cx<-x[i]
   cy<-y[i]
   edist<-sqrt((cx-x)^2+(cy-y)^2)
   temp<-order(edist)
   nbr<-temp[1:k]
   znbr<-z[nbr]
   medz[i]<-median(znbr)
 }
 return(medz)
}


qreq2d <- function(grid,z,zinit,p,cutoff) {
   znew <- zinit
   pznew <- z
   err<-1
   niter<-0
   while ((err>0.005)&(niter<30)) {
       zold <- znew
       pzold <- pznew
       resid <- pzold-zold
       pznew <- zold + dm.rho(resid,p,cutoff)/0.1
   tmp <- as.vector(t(pznew))
       temp <- Tps(grid, tmp)
       znew <- matrix(temp$fitted.values,ncol=ncol(pznew),byrow=T)
       #err<-mean(abs(znew-zold)/pmax(abs(znew),0.01))
   n <- nrow(znew)
   m <- ncol(znew)
   err <- sum((zold-znew)^2)/(n*m)
       niter<-niter+1
#    cat(niter, err, "\n")
   }
   return(znew)
}

qreq2d.wavelet <- function(z,zinit,p,cutoff){
   znew <- zinit
   pznew <- z
   err<-1
   niter<-0
   while ((err>0.05)&(niter<30)) {
   zold <- znew
       pzold <- pznew
       resid <- pzold-zold
       pznew <- zold + 20*dm.rho(resid,p,cutoff)/0.1
   znew <- imwr(threshold(imwd(pznew),type = "hard", policy = "universal"))
   err0<-mean(abs(znew-zold)/pmax(abs(znew),0.01))
   err1 <- mean(abs(zold-znew))
   err2 <- mean((zold-znew)^2)
   err <- min(err0,err1,err2)
       niter<-niter+1
   #    cat(niter, err, "\n")
   }
   return(znew)
}


