#####################################
# DATA PROCESSING #
#####################################

# Read data source
data <- read.csv('NAFZ_4.csv',sep=',',dec=',')

# Take the shocks within longitude 26-40 and latitude 36-42
data <- data[data$longitude<=400,]

# Clean data format
source('codes_cleaning.R')
data <- clean.shocks(data)

# Label main/after/fore-shocks
source('codes_labeling.R')
data <- label.shocks(data)

# Extract mainshock data with largest aftershock
source('codes_extract.R')
ms_data <- extract(data)


#####################################
# DESCRIPTIVE ANALYSIS #
#####################################

set.seed(800)


pdf('plot/ts_shocks.pdf',7,5)
plot(data$date,data$magnitude,pch=16,cex=.6)
dev.off()

pdf('plot/ts_shocks_lab.pdf',7,5)
plot(data$date,data$magnitude,pch=16,cex=.6)
points(ms_data$date,ms_data$magnitude,pch=16,cex=.6,col='white')
points(ms_data$date,ms_data$magnitude,pch=4,cex=.7,col='red')
dev.off()

## Truncate after 1965

tp2 <- ISOdate(1965,01,01,00,00,00,tz="GMT")
data <- data[data$date >tp2,]
ms_data <- ms_data[ms_data$date>tp2,]

pdf('plot/ts_shocks2.pdf',7,5)
plot(data$date,data$magnitude,pch=16,cex=.6,xlab='date',ylab='magnitude')
dev.off()

pdf('plot/ts_shocks_lab2.pdf',7,5)
plot(data$date,data$magnitude,pch=16,cex=.6,xlab='date',ylab='magnitude')
points(ms_data$date,ms_data$magnitude,pch=16,cex=.6,col='white')
points(ms_data$date,ms_data$magnitude,pch=4,cex=.7,col='red')
legend('topright',pch=c(4),col=c('red'),c('main'))
dev.off()


###################################################
## inspect time stationarity ##

attach(ms_data)

par(mfrow=c(1,1))

pdf('plot/ts_ms.pdf',7,3.5)
plot(date,magnitude,pch=16,cex=.7,ylab='mainshock')
dev.off()
pdf('plot/ts_mas.pdf',7,3.5)
plot(date,as_max,pch=16,cex=.7,ylab='max aftershock')
dev.off()

#####################################
# PARAMETRIC ANALYSIS #
#####################################

## fit an exponential distribution to the mainshocks ##
ms <- ms_data$magnitude
ms.jitter <- ms+runif(length(ms),-.05,.05)
# hist(ms)
require(vcd)
require(MASS)
# Use 5 as location end point
fit <- fitdistr(ms-4.95, "exponential")
pdf('plot/fit_ms.pdf',7,5)
hist(ms.jitter,freq=F,breaks=30,main=NULL,xlab='mainshock, jittered')
curve(dexp(x-4.95,fit$estimate),add=T)
legend('topright',c('fitted density'),lty=1)
dev.off()
# Test for goodness of fit
ks.test(ms.jitter-4.95, "pexp", fit$estimate)

# QQ plot
qqplot(ms.jitter,rexp(length(ms.jitter),fit$estimate)+4.95,xlab='mainshock',
       ylab='exponential')
abline(0,1,col='red',lty=2)

############################################################
## fit a Gompertz distribution to Mainshock - MaxAftershock ##
# The aftershock non-truncated obs
as <- ms_data$as_max
as1 <- as[!is.na(as)]
ms1 <- ms[!is.na(as)]
# The aftershock truncated obs
ms2 <- ms[is.na(as)]
# Histogram of the difference mainshock - aftershock, jittered.  
# Note that this might not be the desired Gompertz distribution as those with smaller aftershocks are truncated.
hist(ms1-as1+runif(length(ms1),-.05,.05),breaks=20, main='mainshock - aftershock',xlab='magnitude')

# fitting part I -- the non-truncated obs
x <- ms1-as1
n.x <- length(x)
x2 <- ms2-3.95
# b: scale >0
# eta: shape >0
ll_gom <- function(par.vec){
  b <- par.vec[1]
  eta <- par.vec[2]
  # n.x*log(b) + n.x*log(eta) + n.x*eta + sapply(b,function(b)sum(x)*b) - sapply(eta,function(eta){
  #   eta*sapply(b,function(b)sum(exp(b*x)))
  # })
  # part 1 non-truncated    # part 2 truncated
  -n.x*log(b) - n.x*log(eta) - n.x*eta - sum(x)*b + eta*sum(exp(b*x)) + eta*sum(exp(b*x2)-1) 
  
}
b <- optim(c(1,.4), ll_gom)$par[1]
eta <- optim(c(1,.4), ll_gom)$par[2]

############################################################
## Plot the fit ##

# plot the histogram of truncated MS-AS with the fitted gompertz
dgom <- function(x){
  b*eta*exp(eta)*exp(b*x)*exp(-eta*exp(b*x))
}
# hist(ms1-as1+runif(length(ms1),-.05,.05),freq=F,breaks=20, main='mainshock - aftershock', xlab='magnitude',
#      xlim=c(0,4))
# curve(dgom,xlim=c(0,4),add=T)
# testing the goodness of fit
cgom <- function(x){
  1- exp(-eta*(exp(b*x)-1))
}
# ks.test(ms1-as1+runif(length(ms1),0,.1),cgom)

# randomly generating from this Gompertz distribution with truncation
rgom <- function(t){
  r <- runif(100,0,1)
  x <- 1/b*log(1-log(1-r)/eta)
  x[x>t][1]
}

# Fill in the missing aftershocks
as2 <- ms2 - sapply(ms2-3.95,rgom)

# plot the existing aftershocks
pdf('plot/fill_missing.pdf',7,5)
plot(ms1+runif(length(ms1),-.05,.05),as1+runif(length(ms1),-.05,.05),xlab='mainshock, jittered',ylab='aftershock, jittered',pch=16,cex=.5,
     ylim=c(0,6.5))
# plot(ms1,as1,xlab='mainshock, jittered',ylab='aftershock, jittered',pch=16,cex=.5,
# ylim=c(0,6.5))
points(ms2+runif(length(ms2),-.05,.05),as2,pch=16,cex=.5,col='red')
# points(ms2,as2,pch=16,cex=.5,col='red')
abline(0,1)
abline(-1.1,1,lty=2)
legend('bottomright',lty=c(NA,NA,1,2),pch=c(16,16,NA,NA),col=c('black','red','black','black'),
       c('observed','simulated censored','y=x','y=x-1.1'))
dev.off()

# Put both axes at the same scale
pdf('plot/fill_missing.pdf',6,6)
plot(ms1+runif(length(ms1),-.05,.05),as1+runif(length(ms1),-.05,.05),xlab='mainshock, jittered',ylab='aftershock, jittered',pch=16,cex=.5,
     ylim=c(2.5,6),xlim=c(4.5,8))
# plot(ms1,as1,xlab='mainshock, jittered',ylab='aftershock, jittered',pch=16,cex=.5,
# ylim=c(0,6.5))
points(ms2+runif(length(ms2),-.05,.05),as2,pch=16,cex=.5,col='red')
# points(ms2,as2,pch=16,cex=.5,col='red')
abline(0,1)
abline(-1.1,1,lty=2)
legend('bottomright',lty=c(NA,NA,1,2),pch=c(16,16,NA,NA),col=c('black','red','black','black'),
       c('observed','simulated censored','y=x','y=x-1.1'))
dev.off()

jit.noise <- runif(length(ms1),-.05,.05)
pdf('plot/fit_diff.pdf',7,5)
hist(c(ms1-as1+jit.noise,ms2-as2),freq=F,breaks=20,
     main=NULL,xlab='mainshock - max aftershock, jittered')
curve(dgom,xlim=c(0,6),add=T)
legend('topright',c('fitted density'),lty=1)
dev.off()
ks.test(c(ms1-as1+jit.noise,ms2-as2),cgom)


# QQ plot
r <- runif(length(ms.jitter),0,1)
sim.noise <- 1/b*log(1-log(1-r)/eta)
qqplot(c(ms1-as1+jit.noise,ms2-as2),sim.noise,xlab='Z',ylab='Gompertz distribution')
abline(0,1,col='red',lty=2)

######################################################
## Calculate exceedance probability ##
# MS>x, AS>y
f.ms <- function(x){
  dexp(x-4.95,fit$estimate)
}
f.diff <- function(z){
  b*eta*exp(b*z)*exp(eta)*exp(-eta*exp(b*z))
}
lambda <- fit$estimate 
tmp_fn2 <- function(m){
  exp(-lambda*(m-4.95))
}
tmp_f2 <- function(m){
  sapply(m,tmp_fn2)
}
cond_prob <- integrate(tmp_f2,4.95,Inf)$value
ex_prob <- function(x,y){
  # tmp_fn <- function(z){
  #   exp(-lambda*(y+z))*b*eta*exp(b*z)*exp(eta)*exp(-eta*exp(b*z))
  # }
  tmp_f <- function(z){
    sapply(z,dgom)
  }
  tmp_fn2 <- function(m){
    integrate(tmp_f,0,m-y)$value * exp(-lambda*(m-4.95))
  }
  tmp_f2 <- function(m){
    sapply(m,tmp_fn2)
  }
  integrate(tmp_f2,x,Inf)$value/cond_prob
}

ex_prob(6.5,5.8)
ex_prob(7.1,5.2)
ex_prob(10,8)
ex_prob(10,6)
ex_prob(10.85,6.24)
ex_prob(4.95,0)

#############################################

## TO BE REVISED ##

## Plot level curves ##

source('juan.r')

pdf('plot/level.pdf',8,6)
plot(c(4,11),c(4,11), xlim=c(5,12),type='n',xlab='mainshock',ylab='max aftershock',main='Level Curves')
points(XY,pch=16,cex=.7);
q.list <- c(0.001,0.0005,0.0001,0.00005,0.00001,0.000005,0.000001)
# q.list <- seq(0.001,0.0001,by=-0.0001)
for (i in 1:length(q.list)){
  q <- q.list[i]
  x.thres1 <- uniroot(function(x){ex_prob(x,4)-q},interval=c(5,15))$root
  x.thres2 <- uniroot(function(x){ex_prob(x,x)-q},interval=c(6,15))$root
  print(x.thres2)
  lp <- function(x){
    uniroot(function(y){ex_prob(x,y)-q},interval=c(4.01,x))$root
  }
  lp.fn <- function(x){sapply(x,lp)}
  curve(lp.fn,xlim=c(x.thres2+0.001,x.thres1-0.001),col='red',lty=2,lwd=1.5,add=T)
  # lines(c(x.thres1-0.0001,x.thres1),c(lp.fn(x.thres1-0.0001),4))
}
test=lcb(p0=c(1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6))
abline(0,1,lty=2);
abline(-1.1,1,lty=2)
legend('bottomright',lty=c(1,2),lwd=c(1,1.5),col=c('black','red'),
       c('EVT','parametric'))
dev.off()

####################################
## Large earthquake return period ##

attach(ms_data)
disaster <- ms_data[order(magnitude,decreasing=T)[1:10],]
disaster
attach(disaster)
prob_pm <- c()
for (i in 1:10){
  prob_pm <- c(prob_pm,ex_prob(disaster$magnitude[i],disaster$as_max[i]))
}
prob_evt <- c()
for (i in 1:10){
  prob_evt <- c(prob_evt,tpb(disaster$magnitude[i],disaster$as_max[i]))
}
location <- c('Izmit','Gediz','Duzce','Mudurnu','Erzincan',
              'Afyon','Alasehir','Bartin','Dinar','Mugla Province')
dis_table <- data.frame(date=as.character(date),type=type,ms=magnitude,mas=as_max,
                        prob_pm=prob_pm,prob_evt=prob_evt,location=location)
library(xtable)
xtable(dis_table,digits=5)


