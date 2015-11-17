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
# See also:    twopop.R, data.R
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
source("loss.R")

## descriptive plot
# plot kt of 36 countries including China
plot (kt.Sweden.female, type = "l", ylim = c (-250,150), xlab = "Time", ylab ="kt")
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
  t.temp = eval (parse (text = paste (names[i])))
  t = t.temp$year
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
plot (Sweden$year, sm.kt.Sweden.female, type = "l", ylim = c (-250,150), xlab = "Time", ylab = "kt")
for (i in 1: (loop1 - 2))
{
  lines (eval (parse (text = paste (names[i])))$year,
        eval (parse (text = paste ("sm.kt.", names[i], ".female", sep = ""))), col = i)
}
lines (years.mort, sm.kt.China.female, col = "black", lwd = 3)

#### starting values of thetas

### set up the initial reference curve based on 17 countries Austrila, Austria, Bulgaria,
### Canada, Denmark, Finland, France, Iceland, Italy, Japan, Netherland, Norway,
### Spain, Sweden, Switzerland, UK and USA

names17 = c ("Australia","Austria","Bulgaria","Canada",
        "Denmark","Finland","France","Iceland",
        "Italy","Japan","Netherlands","Norway",
        "Spain","Switzerland",
        "UnitedKingdom","USA","Sweden")
merge1 = sm.kt.Australia.female
for (i in 1:16)
{
  nam14 = paste ("merge", i+1, sep = "")
  temp3 = merge.zoo (eval (parse (text = paste("merge", i, sep = ""))), eval (parse (text = paste ("sm.kt.", names17[i+1], ".female", sep = ""))))
  assign (nam14, temp3)
}
reference0 = rowMeans (merge17, na.rm = TRUE)
reference = ts (reference, start = 1751, frequency = 1)

## plot the reference curve among all 36 smoothed curves
plot (Sweden$year, sm.kt.Sweden.female, type = "l", ylim = c (-250,150), xlab = "Time", ylab = "kt")
for (i in 1: (loop1 - 2))
{
  lines (eval (parse (text = paste (names[i])))$year,
        eval (parse (text = paste("sm.kt.", names[i], ".female", sep = ""))), col = i)
}
lines (years.mort, sm.kt.China.female, col = "black", lwd = 3)
lines (reference, lwd = 4, col = "red")

### find the optimal initial theta based on the reference curve

loss = function (theta, t, kt, t.reference, kt.reference) {
  theta1 = theta[1]
  theta2 = theta[2]
  theta3 = theta[3]
  theta4 = theta[4]
  dref = data.frame (kt.reference, t.reference)
  sm.t = (t - theta2) /theta3 # time adjustment
  sm = locpol (kt.reference~t.reference, dref, kernel = EpaK, xeval = sm.t) # time-adjusted kt based on smoothed reference curve
  mu = theta1 * sm$lpFit[,2] + theta4 # modelled new kt
  mse = mean ((kt-mu)^2) # mse of new kt and the smoothed one
  return (mse)
}

t.reference = 1751:2014
kt.reference = reference
for (i in 1: (loop1 - 1))
{
  # nonlinear optimization
  #if (max (eval (parse (text = paste (names[i])))$year) <= 2011)
  #{theta0 = c (1,0,1,0)}
  #else {if (max (eval (parse (text = paste (names[i])))$year) == 2012)
    #{theta0 = c(1,1,1,0)}
  #else {if (max (eval (parse (text = paste (names[i])))$year) == 2013)
    #{theta0 = c (1,2,1,0)}
  #else {theta0 = c(1,3,1,0)}}}
  theta0 = c (1,0,1,0)
  out = optim (theta0, loss, gr = NULL, eval (parse (text = paste (names[i])))$year,
            eval (parse (text = paste ("sm.kt.", names[i], ".female", sep = ""))),
            t.reference, kt.reference, control = list (maxit = 1000))
  nam7 = paste ("theta0", names[i], "female", sep = ".")
  assign (nam7, out$par)
}
## initial thetas from China
out.China = optim (theta0, loss, gr = NULL, years.mort, sm.kt.China.female,
          t.reference, kt.reference, control = list (maxit = 1000))
theta0.China.female = out.China$par

## test of theta (need set up criteron for next loop)
for (i in 1: loop1)
{
  nam15 = paste ("error.theta", names[i], "female", sep = ".")
  temp4 = mean ((eval (parse (text = paste ("theta0", names[i], "female", sep = ".")))[1] - theta0[1])^2,
             (eval (parse (text = paste ("theta0", names[i], "female", sep = ".")))[2] - theta0[2])^2,
          (eval (parse (text = paste ("theta0", names[i], "female", sep = ".")))[3] - theta0[3])^2,
          (eval (parse (text = paste ("theta0", names[i], "female", sep = ".")))[4] - theta0[4])^2)
  assign (nam15, temp4)
}

## test (the shifted kt's vs previous one) 
loss = function (theta, t, kt, t.reference, kt.reference) {
  theta1 = theta[1]
  theta2 = theta[2]
  theta3 = theta[3]
  theta4 = theta[4]
  dref = data.frame (kt.reference, t.reference)
  sm.t = (t - theta2) /theta3 # time adjustment
  sm = locpol (kt.reference~t.reference, dref, kernel = EpaK, xeval = sm.t) # time-adjusted kt based on smoothed reference
  mu = theta1 * sm$lpFit[,2] + theta4 # modelled new kt
  return (mu)
}

for (i in 1: (loop1 - 1))
{
  nam13 = paste ("test", names[i], sep = ".")
  assign (nam13, loss (eval (parse (text = paste ("theta0", names[i], "female", sep = "."))), eval (parse (text = paste (names[i])))$year,
                    eval (parse (text = paste ("sm.kt.", names[i], ".female", sep = ""))),
                    t.reference, kt.reference))
}
test.China = loss (theta0.China.female, years.mort, sm.kt.China.female, t.reference, kt.reference )
## plot shifted kt of 36 countries including China
plot (reference, type = "l", ylim = c (-250,150), xlab = "Time", ylab = "kt", lwd = 4, col = "red")
for (i in 1: (loop1 -1)) 
{
  lines (eval (parse (text = paste (names[i])))$year,
        eval (parse (text = paste ("test", names[i], sep = "."))), col = i)
}
lines (years.mort, test.China, col = "black", lwd = 3)

## error between shifted curve and previous one (need set up criteron for next loop)
for (i in 1: loop1)
{
  nam16 = paste ("error.curve", names[i], "female", sep = ".")
  temp5 = mean ((eval (parse (text = paste ("test", names[i], sep = "."))) - eval (parse (text = paste ("sm.kt.", names[i], ".female", sep = ""))))^2)
  assign (nam16, temp5)
}









# standardize theta
theta.matrix = matrix (rep (0, 144), 36, 4)
theta.matrix[1,] = theta0.Australia.female
for (i in 2: loop1)
{
  theta.matrix[i,] = theta.matrix[i-1,] + eval (parse (text = paste ("theta0", names[i], "female", sep = ".")))
  }
theta.temp = theta.matrix[36,]
for (i in 1: loop1)
{
  nam8 = paste ("theta0", names[i], "female1", sep = ".")
  assign (nam8, eval (parse (text = paste ("theta0", names[i], "female", sep = ".")))[1] / theta.temp[1])
  nam9 = paste ("theta0", names[i], "female3", sep = ".")
  assign (nam9, eval (parse (text = paste ("theta0", names[i], "female", sep = ".")))[3] / theta.temp[3])
  nam10 = paste ("theta0", names[i], "female2", sep = ".")
  assign (nam10, eval (parse (text = paste ("theta0", names[i], "female", sep = ".")))[2] - theta.temp[2]/ 32)
  nam11 = paste ("theta0", names[i], "female4", sep = ".")
  assign (nam11, eval (parse (text = paste ("theta0", names[i], "female", sep = ".")))[4] - theta.temp[4]/ 32)
}
#construct initial common trend
g = function (theta2, theta3, t, kt) {
  d = data.frame (kt, t)
  sm.t = theta3 * t + theta2 # time adjustment
  sm = locpol (kt~t, d, kernel = EpaK, xeval = sm.t) 
  mu = sm$lpFit[, 2]
  return (mu)
}

for (i in 1: (loop1 -1))
{
  nam12 = paste ("g", i, sep = "")
  assign (nam12, g (eval (parse (text = paste ("theta0", names[i], "female", sep = ".")))[2],
                 eval (parse (text = paste ("theta0", names[i], "female", sep = ".")))[3], eval (parse (text = paste (names[i])))$year,
                 eval (parse (text = paste ("kt.", names[i], ".female", sep = "")))))
}
 
