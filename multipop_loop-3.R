# ------------------------------------------------------------------------------
# Project:     Mortality Model for Multip-populations: A Semiparametric 
#              Comparison Approach
# ------------------------------------------------------------------------------
# Quantlet:    multipop.R
# ------------------------------------------------------------------------------
# Description: Estimate and forecast mortality rates based on a semi-parametric 
#              approach, which applies parametric modelling for multiple 
#              nonparametric curves with the shape-related nonlinear variation.
# ------------------------------------------------------------------------------
# Keywords:    nonparametric smoothing, parametric modeling, common trend,
#              mortality, Lee-Carter method, multi-populations
# ------------------------------------------------------------------------------
# See also:    twopop.R, data.R, optimization.R, referencecurve.R, 
#              normalization.R
# ------------------------------------------------------------------------------
# Author:      Lei Fang
# ------------------------------------------------------------------------------

rm(list=ls(all=TRUE))
# set the work directory
setwd("/Users/lei.fang/Desktop/multi-populations model/semipop R")

# install packages
library (demography)
library (locpol)
library (rgl)
library(KernSmooth)
library(sm)
library(boot)

par (mar = c (5, 5, 2, 2), cex.axis = 1.5, cex.lab = 2)

#source("data.R")
load("~/Desktop/multi-populations model/semipop R/data.RData")
source("optimization-3.R")
source("normalization-3.R")
source("referencecurve.R")

# read multi-pop female mortality of 35 countries from Human Mortality Database
shortnames.all = c ("AUS","AUT","BLR","BGR","CAN","CHL","CZE","DNK","EST","FIN","FRATNP",
                    "DEUTNP","HUN","ISL","IRL","ISR","ITA","JPN","LVA","LTU","LUX","NLD","NZL_NP",
                    "NOR","POL","PRT","RUS","SVK","SVN","ESP","CHE","TWN","GBR_NP","USA","SWE")
names.all = c ("Australia","Austria","Belarus","Bulgaria","Canada","Chile","CzechRepublic",
               "Denmark","Estonia","Finland","France","Germany","Hungary","Iceland","Ireland","Israel",
               "Italy","Japan","Latvia","Lithuania","Luxembourg","Netherlands","NewZealand","Norway",
               "Poland","Portugal","Russia","Slovakia","Slovenia","Spain","Switzerland",
               "Taiwan","UnitedKingdom","USA","Sweden","China")

loop.all = length(names.all)

# set bandwidth
bw.default = 2.5

## descriptive plot
# plot kt of 36 countries including China
plot (kt.Sweden.female, type = "l", ylim = c (-250,150), xlab = "Time", ylab ="kt", 
      main = "Original Kt (36 countries)")
for(i in 1: (loop.all - 2))
{
  lines (eval (parse (text = paste ("kt.", names.all[i], ".female", sep = ""))), col = i)
}
lines (kt.China.female, col = "black", lwd = 3)

##### common trend

#### initial setting
### nonparametric smoothing 36 countries including China
for(i in 1: (loop.all - 1))
{
  kt = eval (parse (text = paste ("kt.", names.all[i], ".female", sep = "")))
  t = eval (parse (text = paste (names.all[i])))$year
  d = data.frame(kt,t)
  
  # way 1 - locpol: locpol
  #sm = locpol (kt~t, d, kernel = EpaK, xeval = t, bw = bw.default) # smooth kt
  #nam6 = paste ("sm.kt", names.all[i], "female", sep = ".")
  #temp6 = ts (sm$lpFit[,2], start = t[1], frequency = 1)
  #assign (nam6, temp6)
  
  # way 2 - locpol with Cross Validation bandwidth
  #nam7 = paste ("cvBwSel", names.all[i], "female", sep = ".")
  #temp7 = regCVBwSelC(t,kt, deg = 1, kernel = EpaK, interval = c(0, 10))
  #assign (nam7, temp7)
  #sm = locpol (kt~t, d, kernel = EpaK, xeval = t, bw = temp7) # smooth kt
  #nam6 = paste ("sm.kt", names.all[i], "female", sep = ".")
  #temp6 = ts (sm$lpFit[,2], start = t[1], frequency = 1)
  #assign (nam6, temp6)
  
  # way 3 - sm: sm.regression with optimal smoothing parameter
  h.optimal1 = h.select (t, kt)
  sm = sm.regression(t, kt, h = h.optimal1, eval.points = t, model = "none", poly.index = 1, display="none")
  nam6 = paste ("sm.kt", names.all[i], "female", sep = ".")
  temp6 = ts (sm$estimate, start = t[1], frequency = 1)
  assign (nam6, temp6)
  
  # way 4 - KernSmooth: locpoly
  #sm = locpoly (t, kt, kernel = EpaK, bandwidth = bw.default, gridsize = length(t), range.x = range(t))
  #nam6 = paste ("sm.kt", names.all[i], "female", sep = ".")
  #temp6 = ts (sm$y, start = t[1], frequency = 1)
  #assign (nam6, temp6)
  
  # way 5 - stats: ksmooth
#   sm = ksmooth (t, kt, kernel = "normal", bandwidth = bw.default, n.points  = length(t), range.x = range(t))
#   nam6 = paste ("sm.kt", names.all[i], "female", sep = ".")
#   temp6 = ts (sm$y, start = t[1], frequency = 1)
#   assign (nam6, temp6)
  
  
}

## smooth China female data
d = data.frame (kt.China.female, years.mort)

# way 1 - locpol: locpol
#sm = locpol (kt.China.female~years.mort, d, kernel = EpaK, xeval = years.mort, bw = bw.default)
#sm.kt.China.female = ts (sm$lpFit[,2], start = 1994, frequency = 1)

# way 2 - locpol with Cross Validation bandwidth
#cvBwSel.China.female = regCVBwSelC(years.mort,kt.China.female, deg = 1, kernel = EpaK, interval = c(0, 10))
#sm = locpol (kt.China.female~years.mort, d, kernel = EpaK, xeval = years.mort, bw = cvBwSel.China.female)
#sm.kt.China.female = ts (sm$lpFit[,2], start = 1994, frequency = 1)

# way 3 - sm: sm.regression with optimal smoothing parameter
h.optimal2 = h.select (years.mort, kt.China.female)
sm = sm.regression(years.mort, kt.China.female, h = h.optimal2, eval.points = years.mort, model = "none", poly.index = 1, display="none")
sm.kt.China.female = ts (sm$estimate, start = 1994, frequency = 1)

# way 4 - KernSmooth: locpoly
#sm = locpoly (years.mort, kt.China.female, kernel = EpaK, bandwidth = bw.default, gridsize = length(years.mort), range.x = range(years.mort))
#sm.kt.China.female = ts (sm$y, start = 1994, frequency = 1)

# way 5 - stats: ksmooth
# sm = ksmooth (years.mort, kt.China.female, kernel = "normal", bandwidth = bw.default, n.points  = length(years.mort), range.x = range(years.mort)) 
# sm.kt.China.female = ts (sm$y, start = 1994, frequency = 1)


## plot smoothed kt of 36 countries including China
plot (sm.kt.Sweden.female, type = "l", ylim = c (-250,150), xlab = "Time", ylab = "kt", 
      main = "Smoothed Kt (36 countries)")
for (i in 1: (loop.all - 2))
{
  lines (eval (parse (text = paste ("sm.kt.", names.all[i], ".female", sep = ""))), col = i)
}
lines (sm.kt.China.female, col = "black", lwd = 3)

#### intial values of thetas
### remove Russia, Lithuavia, Latvia, Estonia and Belarus
### set up the initial reference curve based on 31 countries 

names.31 = c ("Australia","Austria","Bulgaria","Canada","Chile","CzechRepublic",
           "Denmark","Finland","France","Germany","Hungary","Iceland","Ireland","Israel",
           "Italy","Japan","Luxembourg","Netherlands","NewZealand","Norway",
           "Poland","Portugal","Slovakia","Slovenia","Spain","Switzerland",
           "Taiwan","UnitedKingdom","USA","Sweden","China")
loop.31 = length (names.31)
merge1 = sm.kt.Australia.female
for (i in 1: (loop.31 -1))
{
  nam7 = paste ("merge", i+1, sep = "")
  temp3 = merge.zoo (eval (parse (text = paste("merge", i, sep = ""))), 
                     eval (parse (text = paste ("sm.kt.", names.31[i+1], ".female", sep = ""))))
  assign (nam7, temp3)
}
reference0.temp = rowMeans (merge31, na.rm = TRUE)

reference0 = ts (reference0.temp, start = Sweden$year[1], frequency = 1)

# reference0.temp.nonsmooth = ts (reference0.temp, start = Sweden$year[1], frequency = 1)
# t.reference0.temp = time(reference0.temp.nonsmooth)
# 
# ## smooth the initial reference curve
# 
# # way 3 - sm: sm.regression
# reference0.temp.smooth =  sm.regression(t.reference0.temp, reference0.temp.nonsmooth, eval.points = t.reference0.temp, model = "none", poly.index = 1, display="none")
# reference0 = ts (reference0.temp.smooth$estimate, start = t.reference0.temp[1], frequency = 1)
# 
# # way 5 - stats: ksmooth
# # reference0.temp.smooth = ksmooth(t.reference0.temp, reference0.temp.nonsmooth, kernel = "normal", bandwidth = bw.default, n.points  = length(t.reference0.temp), range.x = range(t.reference0.temp))
# # reference0 = ts (reference0.temp.smooth$y, start = t.reference0.temp[1], frequency = 1)


## plot the reference curve among all 31 smoothed curves
plot (sm.kt.Sweden.female, type = "l", ylim = c (-250,150), xlab = "Time", ylab = "kt", 
      main = "Reference Curve vs. Smoothed Kt (31 countries)")
for (i in 1: loop.31)
{
  lines (eval (parse (text = paste("sm.kt.", names.31[i], ".female", sep = ""))), col = i)
}
lines (reference0, lwd = 4, col = "red")


##### begin loop 
theta0 = matrix(rep(c (1,0,1),loop.31),loop.31,3,byrow = TRUE)
iteration = 5
for (j in 1 : iteration)
{
### find the optimal initial theta based on the reference curve
theta = eval (parse (text = paste ("theta", j-1, sep = "")))
kt.reference = eval (parse (text = paste ("reference", j-1, sep = "")))
results=matrix(NA,loop.31,5) # c(theta1,theta2,theta3, loss, convergence)
for ( i in 1:loop.31) {
  kt = eval (parse (text = paste ("sm.kt.", names.31[i], ".female", sep = "")))
  tt = time(kt)
  ltt = length(tt)
  temp4 = optimization (theta[i,], kt, kt.reference)
  results[i,] = temp4[[1]]
  nam12 = paste( "shift.kt", names.31[i], j, sep = ".")
  assign (nam12, temp4[[3]])
  nam15 = paste( "mu", names.31[i], j, sep = ".")
  assign (nam15, temp4[[2]])
  pdf(paste("plot_",j,"_",i,".pdf",sep=""))
  par(mar=c(5, 5, 2, 2))
  plot(eval (parse (text = paste ("kt.", names.31[i], ".female", sep = ""))),ylim = c(-250,150), 
       xlim = c(1750,2020),type ="p", xlab = "Time", ylab = "kt",
       main=paste(names.31[i], tt[1], '--', tt[ltt], ':step', j))
  lines(kt, col=2)  # sm.kt
  lines(kt.reference, col=4, lwd=2) # reference
  lines(temp4[[2]], col=5, lty=5) #
  lines(temp4[[3]], col=5, lwd=2) #
  
  abline(v=temp4[[4]], col="gray", lty=5)
  abline(v=temp4[[5]], col="gray", lty=5)
  #legend("topright", legend=c("kt", "sm.kt", "k0", "kt.hat"), col=c(1,2,4,5), lty=1)
  dev.off()
}
nam16 = paste ( "results", j, sep = "")
assign( nam16, results)
nam18 = paste ( "theta", j, sep = "")
assign( nam18, results[,1:3])

## plot shifted kt of 31 countries and the reference curve of previous step
plot (kt.reference, type = "l", ylim = c (-250,150), xlab = "Time", ylab = "kt", lwd = 4, col = "red", 
      main = paste('Reference Curve vs. Shifted Kt (31 countries): step', j-1 ))
for (i in 1: loop.31) 
{
  lines (eval (parse (text = paste ("shift.kt", names.31[i], j, sep = "."))), col = i)
}


# construct optimal.theta.matrix and normalize optimal theta
optimal.theta.matrix = results[,1:3]
standard.optimal.theta.matrix = normalization(optimal.theta.matrix)

# construct new reference curve
for (i in 1: loop.31)
{
  nam13 = paste ("shift.kt", names.31[i], "standard", sep = ".")
  assign(nam13, referencecurve(standard.optimal.theta.matrix[i,2], 
                               standard.optimal.theta.matrix[i,3],
                               eval (parse (text = paste ("sm.kt.", names.31[i], ".female", sep = ""))),
                               eval (parse (text = paste ("shift.kt", names.31[i], j, sep = ".")))))
}

merge1.1 = shift.kt.Australia.standard
for (i in 1: (loop.31 -1))
{
  nam7 = paste ("merge1", i+1, sep = ".")
  temp3 = merge.zoo (eval (parse (text = paste("merge1", i, sep = "."))), 
                     eval (parse (text = paste ("shift.kt", names.31[i+1], "standard", sep = "."))))
  assign (nam7, temp3)
}
reference.temp = rowMeans (merge1.31, na.rm = TRUE)
reference.temp.nonsmooth = ts (reference.temp, start = min(index(merge1.31)), frequency = 1)
t.temp = time(reference.temp.nonsmooth)

## smooth the updated reference curve

# way 3 - sm: sm.regression
temp17 = sm.regression(t.temp, reference.temp.nonsmooth, eval.points = t.temp, model = "none", poly.index = 1, display="none")
nam17 = paste ( "reference", j, sep = "")
assign( nam17, ts (temp17$estimate, start = t.temp[1], frequency = 1))

# way 5 - stats: ksmooth
# temp17 = ksmooth(t.temp, reference.temp.nonsmooth, kernel = "normal", bandwidth = bw.default, n.points  = length(t.temp), range.x = range(t.temp))
# nam17 = paste ( "reference", j, sep = "")
# assign( nam17, ts (temp17$y, start = t.temp[1], frequency = 1))

## plot the updated reference curve among all 31 shifted curves
plot (eval (parse (text = paste ("reference", j, sep = ""))), lwd = 4, 
      col = "red", type = "l", ylim = c (-250,150), xlab = "Time", ylab = "kt", 
      main = paste ('Updated Reference Curve vs. Shifted Kt (31 countries): step', j))
for (i in 1: loop.31)
{
  lines (eval (parse (text = paste("shift.kt", names.31[i], j, sep = "."))), col = i)
  lines (eval (parse (text = paste("sm.kt.", names.31[i], ".female", sep = ""))), col = "grey")
}
lines(eval (parse (text = paste ("reference", j-1, sep = ""))), col = "blue",lwd=3, lty = 5)

}