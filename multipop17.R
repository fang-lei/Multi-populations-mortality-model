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

par (mar = c (5, 5, 2, 2), cex.axis = 1.5, cex.lab = 2)

source("data.R")
source("optimization.R")
source("normalization.R")
source("referencecurve.R")

## descriptive plot
# plot kt of 36 countries including China
plot (kt.Sweden.female, type = "l", ylim = c (-250,150), xlab = "Time", ylab ="kt", 
      main = "Original Kt (36 countries)")
for(i in 1: (loop1 - 2))
{
  lines (eval (parse (text = paste ("kt.", names[i], ".female", sep = ""))), col = i)
}
lines (kt.China.female, col = "black", lwd = 3)

##### common trend

#### initial setting
### nonparametric smoothing 36 countries including China
for(i in 1: (loop1 - 1))
{
  kt = eval (parse (text = paste ("kt.", names[i], ".female", sep = "")))
  t = eval (parse (text = paste (names[i])))$year
  d = data.frame(kt,t)
  sm = locpol (kt~t, d, kernel = EpaK, xeval = t) # smooth kt
  nam6 = paste ("sm.kt", names[i], "female", sep = ".")
  temp6 = ts (sm$lpFit[,2], start = t[1], frequency = 1)
  assign (nam6, temp6)
}

## smooth China female data
d = data.frame (kt.China.female, years.mort)
sm = locpol (kt.China.female~years.mort, d, kernel = EpaK, xeval = years.mort)
sm.kt.China.female = ts (sm$lpFit[,2], start = 1994, frequency = 1)

## plot smoothed kt of 36 countries including China
plot (sm.kt.Sweden.female, type = "l", ylim = c (-250,150), xlab = "Time", ylab = "kt", 
      main = "Smoothed Kt (36 countries)")
for (i in 1: (loop1 - 2))
{
  lines (eval (parse (text = paste ("sm.kt.", names[i], ".female", sep = ""))), col = i)
}
lines (sm.kt.China.female, col = "black", lwd = 3)

##### pre-test with 17 countries
#### intial values of thetas

### set up the initial reference curve based on 17 countries Austrila, Austria, Bulgaria,
### Canada, Denmark, Finland, France, Iceland, Italy, Japan, Netherland, Norway,
### Spain, Sweden, Switzerland, UK and USA

names17 = c ("Australia","Austria","Bulgaria","Canada",
        "Denmark","Finland","France","Iceland",
        "Italy","Japan","Netherlands","Norway",
        "Spain","Switzerland",
        "UnitedKingdom","USA","Sweden")
merge1 = sm.kt.Australia.female
loop2 = length (names17)
for (i in 1: (loop2 -1))
{
  nam7 = paste ("merge", i+1, sep = "")
  temp3 = merge.zoo (eval (parse (text = paste("merge", i, sep = ""))), 
                     eval (parse (text = paste ("sm.kt.", names17[i+1], ".female", sep = ""))))
  assign (nam7, temp3)
}
reference0 = rowMeans (merge17, na.rm = TRUE)
reference = ts (reference0, start = Sweden$year[1], frequency = 1)

## plot the reference curve among all 17 smoothed curves
plot (sm.kt.Sweden.female, type = "l", ylim = c (-250,150), xlab = "Time", ylab = "kt", 
      main = "Reference Curve vs. Smoothed Kt (17 countries)")
for (i in 1: (loop2 - 1))
{
  lines (eval (parse (text = paste("sm.kt.", names[i], ".female", sep = ""))), col = i)
}
lines (reference, lwd = 4, col = "red")

### find the optimal initial theta based on the reference curve
theta0 = c (1,0,1,0)
kt.reference = reference
for ( i in 1:loop2) {
  kt = eval (parse (text = paste ("sm.kt.", names17[i], ".female", sep = "")))
  temp4 = optimization (theta0, kt, kt.reference)
  nam9 = paste( "loss.mse", names17[i], 0, sep = ".")
  assign (nam9, temp4[[1]])
  nam10 = paste( "optimal.theta", names17[i], 0, sep = ".")
  assign (nam10, temp4[[2]])
  nam11 = paste( "error.theta", names17[i], 0, sep = ".")
  assign (nam11, temp4[[3]])
  nam12 = paste( "shift.kt", names17[i], 0, sep = ".")
  assign (nam12, temp4[[4]])
}

## plot shifted kt of 17 countries and the initial reference curve
plot (kt.reference, type = "l", ylim = c (-250,150), xlab = "Time", ylab = "kt", lwd = 4, col = "red", 
      main = "Reference Curve vs. Shifted Kt (17 countries): step 0")
for (i in 1: loop2) 
{
  lines (eval (parse (text = paste ("shift.kt", names17[i], 0, sep = "."))), col = i)
}

# construct optimal.theta.matrix and normalize optimal theta
l = list()
for ( i in 1:loop2) {
  l[[i]] = c (eval (parse (text = paste ("optimal.theta", names17[i], 0, sep = "."))))
}
optimal.theta.matrix = do.call(rbind,l)
normal.theta.matrix = normalization (optimal.theta.matrix)

# updated reference curve / common trend (call referencecurve.R)

for (i in 1: loop2)
{
  nam15 = paste ("g", i, sep = "")
  assign (nam15, referencecurve (normal.theta.matrix[i, 2],
                    normal.theta.matrix[i, 3], 
                    eval (parse (text = paste("sm.kt.", names17[i], ".female", sep = "")))))
}