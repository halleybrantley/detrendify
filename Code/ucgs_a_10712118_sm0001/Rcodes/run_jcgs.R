library(fields)
library(quantreg)
library(stats)
library(wavethresh)
source("function_jcgs.R")
# it is for 1-d simulation (asymmetric noise case) in Section 4.1
n <- 2000
i <- c(1:n)
x <- i/n
ss <- 1
set.seed(167)
noise <- rnorm(n,0,0.07)
ss <- 1
noise<-rgamma(n,3,1)
y <- sin(10*x)+noise
data <- cbind(x,y)
rx <- sort(data[,1])
ry <- data[,2][order(data[,1])]
k <- 100

tau <- 0.10 #(0.5,0.9)
finit <- as.vector(int.fit.fig(k,n,tau,ry))
fit1.two <- aresult(x,(2000-100+1),0.04,hry1(100,2000,x),hry2(100,tau,2000,y)) #Yu's method
fit1.pro2 <- qreq1d.sreg.fig(x,y,fit1.two,tau,0.01) # the proposed method
temp <- rqss(y~qss(x,constraint="N"),tau=tau,lambda=0.1) 
fit1.rqss <- temp$coef[1]+temp$coef[-1] # Koenker's method 
tr1 <- sin(10*x)+qgamma(tau,3,1)

plot(x,tr1,type="l")
lines(x,fit1.two,col=2)
lines(x,fit1.pro2,col=3)
lines(x[-1],fit1.rqss,col=4)

mean((tr1-fit1.pro2)^2)
mean((tr1-fit1.two)^2)
mean((tr1[-1]-fit1.rqss)^2)

# it is for 2-d simulation in Section 5.1
n <- 10
i <- c(1:n)
x <- i/n
y <- x
z <- sin(3*pi*x)%*%t(cos(pi*y))
m <- ncol(z)
noise <- matrix(rnorm(n*m),ncol=m) # for normal noise
z.noise <- z+noise/3

tau <- 0.5
tr <- z+matrix(rep(qnorm(tau),n*m),ncol=m)/3 # true quantile with normal noise

k <- 3
zinit <- d2int.fit.sy(x,y,z.noise,k,tau) # it is for a initial fit
grid <- NULL
for(i in 1:n){
    tmp <- cbind(rep(x[i],length(y)),y)
    grid <- rbind(grid,tmp)
}
zfit <- qreq2d(grid,z.noise,zinit,tau,0.07)

persp(x, y, zfit, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
zlim=c(min(z,tr,zfit),max(z,tr,zfit)),ltheta = 120, ticktype = "detailed")


# it is for John Lennon image in Figure 6
data(lennon)
m <- ncol(lennon)
n <- nrow(lennon)
snr <- 3
sd <- sqrt(var(as.vector(lennon)))
sigma <- sd/snr 
noise <- matrix(rnorm(n*m,mean=0,sd=sigma),ncol=m)
lennon.noise <- lennon+noise
image(lennon.noise,col=gray((32:0)/32))
i <- c(1:n)
x <- i/n
y <- x
tau <- 0.1 #0.5, 0.9
k <- 3
zinit <- d2int.fit.sy(x,y,lennon+noise,k,tau)
cutoff <- 0.02
fit <- qreq2d.wavelet(lennon+noise,zinit,tau,cutoff) # it is the proposed fit according to tau. 
image(fit,col=gray((32:0)/32))

