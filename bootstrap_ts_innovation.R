rm(list=ls(all=TRUE))
# set the work directory
setwd("/Users/lei.fang/Desktop/multi-populations model/semipop R")

options(digits = 4)
# install packages
library (demography)
library (locpol)
library (rgl)
library(KernSmooth)
library(sm)
library(boot)
library(tseries)

# source("data.R") to reduce code-running time 
# from previous procedure (multi_loop-3.R)
load("~/Desktop/multi-populations model/semipop R/sm_3_constraint_loop5.RData")

# Bootstrap to resample from China's mortality time series 
# to construct the CI for shape deviation parameters theta

# Reference time series curve (unchanged). We use for 
# comparing with resampled China's mortality TS
kt.referencet = reference4 

# Original China's mortality TS, from which we are going 
# to take TS bootstrap
boot.data = kt.China.female 

# Initialize theta values for optimazation
theta = theta4[31,] 

### loss function
loss = function (theta, t, kt, t.reference, kt.reference) 
{
  ## assume theta[1]>0, theta[3]>0, if not, take absolute values
  theta1 = abs(theta[1])
  theta2 = theta[2]
  theta3 = abs(theta[3])
  
  if (theta[1]<0 | theta[3]<0) warning("theta1 or theta3 <0, abs value is used")
  
  sm.t = (t.reference - theta2) /theta3 # time adjustment
  
  ## common domain for kt and time-adjusted kt.reference
  tmin = max(min(t), min(sm.t))
  tmax = min(max(t), max(sm.t)) 
  
  # index for common domain
  i0 = which(t >= tmin & t <= tmax) 
  
  if (length(i0)>0)
  {
    t0 = t[i0]  
    
    ## smooth interpolation of shifted kt.reference on common grid    
    dref = data.frame(sm.t=sm.t, kt.reference=kt.reference)
    
    # way 3 - sm: sm.regression with optimal smoothing parameter
    #       h.optimal3 = h.select(sm.t, kt.reference)
    #       sm = sm.regression(sm.t, kt.reference, h = h.optimal3, eval.points = t0, model = "none", poly.index = 1, display="none")
    #       mu = theta1 * sm$estimate
    
    # way 5 - stats: ksmooth
    sm = ksmooth (sm.t, kt.reference, kernel = "normal", bandwidth = bw.default, n.points = length(t0), range.x = range(t0))
    mu = theta1 * sm$y
    
    ## mean squared error at common grid points
    mse = mean ((kt[i0]-mu)^2) # mse of kt and the modelled one
  } 
  else 
  { 
    mse = 1e9 
  }
  
  return (mse)
}

# Construct the "statistic" function which is needed 
# and very important in 'tsboot'
# tsboot (tseries, statistic, R, l = NULL, sim = "model",...)
boot.func = function ( kt) {
  t = time (kt)
  t.reference = time (kt.reference)
  
  ### parameter estimation
  conv = 1
  theta0 = theta
  
  while(conv!=0) 
  { 
    # check convergence
    out = optim (theta0, loss, gr = NULL, t, kt, t.reference, kt.reference, control = list (maxit = 1000)) # optimization
    conv = out$convergence
    theta0 = out$par
  }  
  # add constraint on theta2
  if (out$par[2]>=-200 & out$par[2]<=200) temp.par = out$par else temp.par = theta
  result=c(temp.par, out$value, out$convergence)
}

# simulate time series with innovation

bootdata.fit = auto.arima(boot.data)
sim = 500
plot (boot.data, xlab = "Time", ylab = "Kt", ylim = c(-150,70), col = "blue", lwd = 10)

for ( i in 1:sim) {
  nam5 = paste ("sim.ts", i, sep = ".")
  # simulate time seriese based on bootstrap of innovation/error term in ARIMA fitting model
  temp5 = simulate(bootdata.fit,nsim = length(boot.data), future = FALSE, bootstrap = TRUE)
  assign(nam5, temp5)
  
  # time series test on plot
  lines(temp5, col = "grey") 
}


# compute parameters theta's based on bootstrapped time series
theta.boot = matrix(rep(c (0, 0, 0, 0, 0), sim), sim, 5, byrow = TRUE)

for ( i in 1:sim) {
  theta.boot[i,] = boot.func(eval (parse (text = paste ("sim.ts", i, sep = "."))))
}

# plot test
hist(theta.boot[,1],xlab = "theta 1",main = "Histogram of theta 1")
quantile(theta.boot[,1])

hist(theta.boot[,2][which(theta.boot[,2]>=-100 & theta.boot[,2]<=100)],xlab = "theta 2",main = "Histogram of theta 2",xlim = c(-5,5),breaks=10)
quantile(theta.boot[,2])

hist(theta.boot[,3][which(theta.boot[,3]>=0.8 & theta.boot[,3]<=1.2)],xlab = "theta 3",main = "Histogram of theta 3",xlim = c(0.93,0.99),breaks=10)
quantile(theta.boot[,3])

# forecast with estimated parameters from bootstrap
boot.forecast = function (theta, kt, kt.reference) 
{
  t = time (kt)
  t.reference = time (kt.reference)
  ### check results on graph
  theta1 = abs(theta[1])
  theta2 = theta[2]
  theta3 = abs(theta[3])
  
  # time adjustment
  sm.t = (t.reference - theta2) /theta3 
  
  ## common grid for kt and shifted kt.reference
  tmin = max(min(t), min(sm.t))
  tmax = min(max(t), max(sm.t)) 
  i0 = which(t >= tmin & t <= tmax)
  t0 = t[i0]  
  
  ## shifted curves (kt.hat)
  dref = data.frame(sm.t = sm.t, kt.reference = kt.reference)
  
  # way 3 - sm: sm.regression with optimal smoothing parameter
  h.optimal4 = h.select (sm.t, kt.reference)
  sm = sm.regression(sm.t, kt.reference, h = h.optimal4, eval.points = sm.t, model = "none", poly.index = 1, display="none")
  mu = theta1 * sm$estimate
  
  # way 5 - stats: ksmooth
  #     sm = ksmooth (sm.t, kt.reference, kernel = "normal", bandwidth = bw.default, n.points = length(sm.t), range.x = range(sm.t))
  #     mu = theta1 * sm$y
  
  mu = ts (mu, start = sm.t[1], frequency = theta3)
  
  # way 3 - sm: sm.regression with optimal smoothing parameter
  h.optimal5 = h.select (sm.t, kt.reference)
  sm0 = sm.regression(sm.t, kt.reference, h = h.optimal5, eval.points = t0, model = "none", poly.index = 1, display="none")
  mu0 = theta1 * sm0$estimate
  
  # way 5 - stats: ksmooth
  #     sm0 = ksmooth (sm.t, kt.reference, kernel = "normal", bandwidth = bw.default, n.points = length(t0), range.x = range(t0))
  #     mu0 = theta1 * sm0$y
  mu0 = ts (mu0, start = t0[1], frequency = 1)
  
  return ( list( mu, mu0, tmin, tmax ))
}  

# forecast
plot( forecast(bootdata.fit, h = 40, level = 95), xlab = "Time", ylab = "Kt", xlim = c(1828, 2050), ylim = c(-200, 150))
lines(kt.referencet, col = "blue", lwd = 5)
lines(boot.data, col = "red", lwd = 5)

for ( i in 1:sim) 
{
  temp = boot.forecast(theta.boot[i, 1:3], boot.data, kt.reference)
  nam12 = paste( "shift.kt.boot", i, sep = ".")
  assign (nam12, temp[[2]])
  nam15 = paste( "mu.boot", i, sep = ".")
  assign (nam15, temp[[1]])
  lines(temp[[2]], col = "green") #
  lines(temp[[1]], col = "grey") #
}

# way 1
# prediction interval
for (i in 2010:2050)
{
  nam = matrix(0,1,500)
for ( j in 1:sim) {
  if (max(time(eval (parse (text = paste ("mu.boot", j, sep = "."))))) >= i) nam[j] = window(eval (parse (text = paste ("mu.boot", j, sep = "."))), i)[1] else nam[j] = NA
}
  nam1 = paste("ts", i, sep = ".")
  assign( nam1, nam)
  nam2 = paste("ts", i, 5,sep = ".")
  assign( nam2, quantile(nam, probs = seq(0, 1, 0.05), na.rm = T)[2])
  nam3 = paste("ts", i, 95,sep = ".")
  assign( nam3, quantile(nam, probs = seq(0, 1, 0.05), na.rm = T)[20])
  nam4 = paste("ts", i, 20,sep = ".")
  assign( nam4, quantile(nam, probs = seq(0, 1, 0.05), na.rm = T)[5])
  nam5 = paste("ts", i, 80,sep = ".")
  assign( nam5, quantile(nam, probs = seq(0, 1, 0.05), na.rm = T)[17])
  #quantile(ts.2010, probs = seq(0, 1, 0.05), na.rm = T)
}

for (j in c(5,20,80,95))
{
  ts.temp = ts(0, start = 2010, end = 2040, frequency = 1)
  for (i in 2010:2040)
  {
    ts.temp[i-2009] = eval (parse (text = paste ("ts", i, j , sep = ".")))
  }
  nam6 = paste("ts.forecast", j,sep = ".")
  assign(nam6, ts.temp)
}

# plot
plot( forecast(bootdata.fit, h = 30, level = 95), xlab = "Time", ylab = "Kt", xlim = c(1828, 2050), ylim = c(-200, 150))
lines(kt.referencet, col = "blue", lwd = 5)
lines(boot.data, col = "red", lwd = 5)

lines(ts.forecast.5, col = "yellow", lwd = 5)
lines(ts.forecast.95, col = "yellow", lwd = 5)

lines(ts.forecast.20, col = "yellow", lwd = 5, lty =2)
lines(ts.forecast.80, col = "yellow", lwd = 5, lty =2)





# way 2
# plot 95% confidence interval
quan.ts.2010 = matrix(0,1,500)
for ( i in 1:sim) {
  quan.ts.2010[i] = window(eval (parse (text = paste ("sim.ts", i, sep = "."))), 2010)[1]
}

quan.index.10 = matrix(0, 1, 500)

for ( i in 1:sim) {
  if (quan.ts.2010[i] >= quantile(quan.ts.2010,probs = c(0.45, 0.55))[1] & quan.ts.2010[i] <= quantile(quan.ts.2010,probs = c(0.45, 0.55))[2])
    quan.index.10[i] = i
}

quan.index.80 = matrix(0,1,500)

for ( i in 1:sim) {
  if (quan.ts.2010[i] >= quantile(quan.ts.2010,probs = c(0.1, 0.9))[1] & quan.ts.2010[i] <= quantile(quan.ts.2010,probs = c(0.1, 0.9))[2])
    quan.index.80[i] = i
}

quan.index.90 = matrix(0, 1, 500)

for ( i in 1:sim) {
  if (quan.ts.2010[i] >= quantile(quan.ts.2010,probs = c(0.05, 0.95))[1] & quan.ts.2010[i] <= quantile(quan.ts.2010,probs = c(0.05, 0.95))[2])
    quan.index.90[i] = i
}

#plot(forecast(bootdata.fit,h=40,level = 95), xlab = "Time", ylab = "Kt", xlim = c(1828,2050),ylim = c(-300,150))
plot(kt.referencet, col = "black", lwd = 5, xlab = "Time", ylab = "Kt", xlim = c(1828, 2040), ylim = c(-300, 150))
lines(boot.data, col = "red", lwd =5)
mu.boot.0 = ts(0,start = 2010, frequency = 1, end = 2010)

for ( i in quan.index.90) {
  #lines(eval (parse (text = paste ("shift.kt.boot", i, sep = "."))), col= "green") #
  lines(window(eval (parse (text = paste ("mu.boot", i, sep = "."))),2008), col= "blue") #
}

for ( i in quan.index.80) {
  #lines(eval (parse (text = paste ("shift.kt.boot", i, sep = "."))), col= "green") #
  lines(window(eval (parse (text = paste ("mu.boot", i, sep = "."))),2008), col= "grey") #
}

for ( i in quan.index.10) {
  #lines(eval (parse (text = paste ("shift.kt.boot", i, sep = "."))), col= "green") #
  lines(window(eval (parse (text = paste ("mu.boot", i, sep = "."))),2008), col= "yellow") #
}

lines(kt.referencet, col = "black", lwd = 5)
lines(mu.China.4, col = 5, lty = 2, lwd = 4)
