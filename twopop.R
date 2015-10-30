setwd("/Users/lei.fang/Desktop/multi-populations model/semipop R")

library(demography)
library(locpol)
library(rgl)

par(mar=c(5, 5, 2, 2),cex.axis = 1.5, cex.lab = 2)
# Japan
# mortality data

japan1<-read.demogdata("Mx_1x1.txt","Exposures_1x1.txt",type="mortality",label="Japan")

## adjusted data
# Japan (restrict age to 90)
japanadjust<-extract.ages(japan1,ages=0:90)

# China
# mortality data

## China male mortility data
data2<-read.table('Chinamortalitymale.txt',header=F,sep='')
data2pop<-read.table('Chinamortalitymalepop.txt',header=F,sep='')
ages.mort<-0:90
years.mort<-1994:2010
china.mort.male<-demogdata(data2,data2pop,ages.mort,years.mort,type='mortality',label='China',name='male')

## China female mortality data
data3<-read.table('Chinamortalityfemale.txt',header=F,sep='')
data3pop<-read.table('Chinamortalityfemalepop.txt',header=F,sep='')
china.mort.female<-demogdata(data3,data3pop,ages.mort,years.mort,type='mortality',label='China',name='female')

# China (presmooth)
china.mort.male.adjust<-smooth.demogdata(china.mort.male,b=65,k=30)
china.mort.female.adjust<-smooth.demogdata(china.mort.female,b=65,k=30)

# kt indicator (general)
# Japan
japan.lca.female<-lca(japanadjust,series="female",adjust="dt")
ax.japan.female.mort=japan.lca.female$ax
bx.japan.female.mort=japan.lca.female$bx

japan.lca.male<-lca(japanadjust,series="male",adjust="dt")
ax.japan.male.mort=japan.lca.male$ax
bx.japan.male.mort=japan.lca.male$bx

# China
china.lca.female<-lca(china.mort.female.adjust,series="female",adjust="dt",max.age=90,interpolate=TRUE)
ax.china.female.mort=china.lca.female$ax
bx.china.female.mort=china.lca.female$bx

china.lca.male<-lca(china.mort.male.adjust,series="male",adjust="dt",max.age=90,interpolate=TRUE)
ax.china.male.mort=china.lca.male$ax
bx.china.male.mort=china.lca.male$bx

# datasets
kt.jm=japan.lca.male$kt
kt.jf=japan.lca.female$kt
kt.cm=china.lca.male$kt
kt.cf=china.lca.female$kt
tj=1947:2012
tc=1994:2010

plot(tj,kt.jf,type="l",col="blue",lwd=4,xlab="Time",ylab="Kt")
lines(kt.jm,col="red",lwd=4)
lines(tc,kt.cf,col="blue",lwd=4)
lines(kt.cm,col="red",lwd=4)

# shift movement of smoothed China and Japan kt (female)
dcf=data.frame(kt.cf,tc)
sm.cf<- locpol(kt.cf~tc,dcf,kernel=EpaK,xeval=tc)
sm.kt.cf=sm.cf$lpFit[,2] # smoothed kt of China

djf=data.frame(kt.jf,tj)
sm.jf<- locpol(kt.jf~tj,djf,kernel=EpaK,xeval=tj)
sm.kt.jf=sm.jf$lpFit[,2] # smoothed kt of Japan

plot(tj,sm.kt.jf,type="l",col="blue",lwd=4,xlab="Time",ylab="Kt")
lines(kt.jf,col="red",lwd=4)
lines(tc,sm.kt.cf,lwd=4)
lines(kt.cf,col="darkgreen",lwd=4)
lines(tc-20,sm.kt.cf,lwd=4,lty=2)
lines(tc-25,sm.kt.cf,lwd=4,lty=3)
lines(tc-23,sm.kt.cf,lwd=4,lty=4)

# time delay
plot(tj,sm.kt.jf,type="l",col="blue",lwd=4,xlab="Time",ylab="Kt")
lines(tc,sm.kt.cf,lwd=4)
lines(tc-23,sm.kt.cf,lwd=4,lty=4)
# vertical shift
plot(tj,sm.kt.jf,type="l",col="blue",lwd=4,xlab="Time",ylab="Kt")
lines(tc,sm.kt.cf,lwd=4)
lines(tc,sm.kt.cf-85,lwd=4,lty=4)


# loss function, smoothing and nonlinear optimization (female)
# loss function + smoothing
loss <- function(theta,tc,kt.cf,tj,kt.jf){
  theta1=theta[1]
  theta2=theta[2]
  theta3=theta[3]
  theta4=theta[4]
  dj=data.frame(kt.jf,tj)
  dc=data.frame(kt.cf,tc)
  sm.t=(tc-theta2)/theta3 # time adjustment
  sm.cf<- locpol(kt.cf~tc,dc,kernel=EpaK,xeval=tc) # smooth kt of China
  sm.kt.cf=sm.cf$lpFit[,2] # smoothed kt.cm
  sm <- locpol(kt.jf~tj,dj,kernel=EpaK,xeval=sm.t) # time-adjusted kt.cm based on smoothed Japan
  mu = theta1*sm$lpFit[,2]+theta4 # modelled new kt of China
  mse = mean((sm.kt.cf-mu)^2) # mse of new kt of China and the smoothed one
  return(mse)
}

# nonlinear optimization
theta0=c(1,23,1,0)
out=optim(theta0, loss, gr=NULL,tc, kt.cf,tj,kt.jf,control = list(maxit=1000))
# out.LB=optim(theta0, loss, gr=NULL,control = list(trace = TRUE,REPORT = 500),method = "L-BFGS-B",lower=c(0.7,1,0.7,-100), upper=c(1.5,47,1.5,50),tc, kt.cm,tj,kt.jm,control = list(maxit=1000))

# goodness of fit
loss <- function(theta,tc,kt.cf,tj,kt.jf){
  theta1=theta[1]
  theta2=theta[2]
  theta3=theta[3]
  theta4=theta[4]
  dj=data.frame(kt.jf,tj)
  dc=data.frame(kt.cf,tc)
  sm.t=(tc-theta2)/theta3 # time adjustment
  sm.cf<- locpol(kt.cf~tc,dc,kernel=EpaK,xeval=tc) # smooth kt of China
  sm.kt.cf=sm.cf$lpFit[,2] # smoothed kt.cm
  sm <- locpol(kt.jf~tj,dj,kernel=EpaK,xeval=sm.t) # time-adjusted kt.cm based on smoothed Japan
  mu = theta1*sm$lpFit[,2]+theta4 # modelled new kt of China
  return(c(sm.t,mu))
}
theta.hat=out$par
kt.cf.new=loss(theta.hat,tc,kt.cf,tj,kt.jf)

## pdf(file="jr1.pdf",bg="transparent")
plot(tj,sm.kt.jf,type="l",col="blue",lwd=4,xlab="Time",ylab="Kt")
lines(kt.jf,col="red",lwd=4)
lines(tc,sm.kt.cf,lwd=4)
lines(kt.cf,col="darkgreen",lwd=4)
lines(kt.cf.new[1:17],kt.cf.new[18:34],type="p",lwd=4)
## dev.off()

# plot loss function
loss <- function(theta1,theta2,theta3,theta4,tc,kt.cf,tj,kt.jf){
  dj=data.frame(kt.jf,tj)
  dc=data.frame(kt.cf,tc)
  sm.t=(tc-theta2)/theta3 # time adjustment
  sm.cf<- locpol(kt.cf~tc,dc,kernel=EpaK,xeval=tc) # smooth kt of China
  sm.kt.cf=sm.cf$lpFit[,2] # smoothed kt.cm
  sm <- locpol(kt.jf~tj,dj,kernel=EpaK,xeval=sm.t) # time-adjusted kt.cm based on smoothed Japan
  mu = theta1*sm$lpFit[,2]+theta4 # modelled new kt of China
  mse = mean((sm.kt.cf-mu)^2) # mse of new kt of China and the smoothed one
  return(mse)
}

# plot loss function of theta2 and theta4
vtheta2=seq(0,50,len=100)
vtheta4=seq(-100,50,len=100)
npar=length(vtheta2)
res=matrix(rep(0,10000),100,100)
for (i in 1:npar){
  for (l in 1:npar){
    res[i,l]=loss(1,vtheta2[i],1,vtheta4[l],tc,kt.cf,tj,kt.jf)}
  res[i,l]=loss(1,vtheta2[i],1,vtheta4[l],tc,kt.cf,tj,kt.jf)
}

persp3d(vtheta2,vtheta4,res,col = "blue",xlab = "theta2", ylab = "theta4", zlab = "Loss")

resnew=res*(res<50)
persp3d(vtheta2,vtheta4,resnew,col = "blue",xlab = "theta2", ylab = "theta4", zlab = "Loss")

# contour plot
contour(vtheta2,vtheta4,resnew,col=rainbow(20),nlevels=70,lwd=1)
contour(vtheta2,vtheta4,res,col=rainbow(20),nlevels=150,lwd=2)


# plot loss function 2 ---detailed of theta2 and theta4
vtheta22=seq(20,26,len=100)
vtheta42=seq(-5,5,len=100)
npar2=length(vtheta22)
res2=matrix(rep(0,10000),100,100)
for (i in 1:npar2){
  for (l in 1:npar2){
    res2[i,l]=loss(1,vtheta22[i],1,vtheta42[l],tc,kt.cf,tj,kt.jf)}
  res2[i,l]=loss(1,vtheta22[i],1,vtheta42[l],tc,kt.cf,tj,kt.jf)
}

persp3d(vtheta22,vtheta42,res2,col = "blue",xlab = "theta2", ylab = "theta4", zlab = "Loss")

resnew1=res2*(res2<50)
persp3d(vtheta22,vtheta42,resnew1,col = "blue",xlab = "theta2", ylab = "theta4", zlab = "Loss")

# contour plot
contour(vtheta22,vtheta42,resnew1,col=rainbow(20),nlevels=90,lwd=1,xlab="theta2",ylab="theta4")

# plot loss function of theta2
vtheta23=seq(0,50,len=100)
npar3=length(vtheta23)
res3=rep(0,100)
for (i in 1:npar3){
  res3[i]=loss(1,vtheta23[i],1,0,tc,kt.cf,tj,kt.jf)
}

plot(vtheta23,res3,col = "blue",xlab = "theta2", ylab = "Loss",type="l",lwd=4)

# forecast of China kt and mortality rate (female)
kt.c.fore <- function(theta,tc.fore,tj,kt.jf){
  theta1=theta[1]
  theta2=theta[2]
  theta3=theta[3]
  theta4=theta[4]
  dj=data.frame(kt.jf,tj)
  sm.t=(tc.fore-theta2)/theta3 # time adjustment
  sm <- locpol(kt.jf~tj,dj,kernel=EpaK,xeval=sm.t) # time-adjusted kt.cm based on smoothed Japan
  kt.fore = theta1*sm$lpFit[,2]+theta4 # modelled new kt of China
  return(kt.fore)
}

tc.fore=2011:2030
kt.cf.fore=kt.c.fore(theta.hat,tc.fore,tj,kt.jf)
plot(china.mort.female.adjust,ylim=c(-10,0))
yt.cf.fore=matrix(rep(0,1820),20,91)
#yt.cf.fore1=ax.china.female.mort+bx.china.female.mort*kt.cf.fore[1]
#plot(0:90,yt.cf.fore1,type="l",lwd=4,col="grey",ylim=c(-10,-2),xlab="Age",ylab="Log death rate")
for (i in 1:20){
  yt.cf.fore[i,]=ax.china.female.mort+bx.china.female.mort*kt.cf.fore[i]
  lines(yt.cf.fore[i,],lwd=4,col="grey")
}

plot(tj,kt.jf,type="l",lwd=4,col="red",xlab="Time",ylab="Kt",xlim=c(1947,2030),ylim=c(-120,150))
lines(tj,sm.kt.jf,lwd=4,col="blue")
lines(tc,kt.cf,lwd=4,col="red")
lines(tc,sm.kt.cf,lwd=4,col="blue")
lines(tc.fore,kt.cf.fore,lwd=4)
