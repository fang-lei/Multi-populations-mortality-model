# set the work directory
setwd("/Users/lei.fang/Desktop/multi-populations model/semipop R")

# install packages
library (demography)
library (locpol)
library (rgl)
library(KernSmooth)
library(sm)
library(boot)

#source("data.R")
load("~/Desktop/multi-populations model/semipop R/sm_3_loop5.RData")

# bootstrap
kt.referencet = reference4
boot.data = kt.China.female
theta=theta4[31,]


boot.func = function ( kt) {
  t = time (kt)
  t.reference = time (kt.reference)
  ### loss function
  loss = function (theta, t, kt, t.reference, kt.reference) {
    ## assume theta[1]>0, theta[3]>0, if not, take absolute values
    theta1 = abs(theta[1])
    theta2 = theta[2]
    theta3 = abs(theta[3])
    
    if (theta[1]<0 | theta[3]<0) warning("theta1 or theta3 <0, abs value is used")
    sm.t = (t.reference - theta2) /theta3 # time adjustment
    ## common domain for kt and time-adjusted kt.reference
    tmin=max(min(t), min(sm.t))
    tmax=min(max(t), max(sm.t)) 
    i0=which(t>=tmin & t<=tmax) # index for common domain
    if (length(i0)>0){
      t0=t[i0]  
      ## smooth interpolation of shifted kt.reference on common grid    
      dref=data.frame(sm.t=sm.t, kt.reference=kt.reference)
      
      # way 3 - sm: sm.regression with optimal smoothing parameter
#       h.optimal3 = h.select(sm.t, kt.reference)
#       sm = sm.regression(sm.t, kt.reference, h = h.optimal3, eval.points = t0, model = "none", poly.index = 1, display="none")
#       mu = theta1 * sm$estimate
      
      # way 5 - stats: ksmooth
              sm = ksmooth (sm.t, kt.reference, kernel = "normal", bandwidth = bw.default, n.points = length(t0), range.x = range(t0))
              mu = theta1 * sm$y
      
      ## mean squared error at common grid points
      mse = mean ((kt[i0]-mu)^2) # mse of kt and the modelled one
    } else { mse=1e9 }
    return (mse)
  }
  ### parameter estimation
  conv=1
  theta0 = theta
  while(conv!=0) { # check convergence
    out = optim (theta0, loss, gr = NULL, t, kt, t.reference, kt.reference, control = list (maxit = 1000)) # optimization
    conv=out$convergence
    theta0=out$par
  }  
  result=c(out$par, out$value, out$convergence)
}
  
boot.1 = tsboot(boot.data, boot.func, R = 50, l = 12, sim = "fixed")
boot.1$t
boot.1$t0
boot.ci(boot.1, type="bca", index=1)
