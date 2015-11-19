# ------------------------------------------------------------------------------
# Project:     Mortality Model for Multip-populations: A Semiparametric 
#              Comparison Approach
# ------------------------------------------------------------------------------
# Quantlet:    referencecurve.R
# ------------------------------------------------------------------------------
# Description: Regenerate reference curve (common trend) based on normalized 
#              optimal theta and original smoothed kt.
# ------------------------------------------------------------------------------
# Keywords:    normalization, mortality, Lee-Carter method, commom trend
# ------------------------------------------------------------------------------
# See also:    twopop.R, multipop.R, data.R, optimization.R
# ------------------------------------------------------------------------------
# Author:      Lei Fang
# ------------------------------------------------------------------------------

referencecurve = function (theta, kt) {
  theta.temp = colSums(theta)
  for (i in 1: loop2)
  {
    nam10 = paste ("normal.theta", names17[i], 1, sep = ".")
    assign (nam10, theta[i, 1] / theta.temp[1])
    nam11 = paste ("normal.theta", names17[i], 3, sep = ".")
    assign (nam11, theta[i, 3] / theta.temp[3])
    nam12 = paste ("normal.theta", names17[i], 2, sep = ".")
    assign (nam12, theta[i, 2] - theta.temp[2]/ loop2)
    nam13 = paste ("normal.theta", names17[i], 4, sep = ".")
    assign (nam13, theta[i, 4] - theta.temp[4]/ loop2)
    nam14 = paste ("normal.theta", names17[i], sep = ".")
    assign (nam14, c (eval (parse (text = paste ("normal.theta", names17[i], 1, sep = "."))),
                      eval (parse (text = paste ("normal.theta", names17[i], 2, sep = "."))),
                      eval (parse (text = paste ("normal.theta", names17[i], 3, sep = "."))),
                      eval (parse (text = paste ("normal.theta", names17[i], 4, sep = ".")))))
  }
  ll = list()
  for ( i in 1:loop2) {
    ll[[i]] = c (eval (parse (text = paste ("normal.theta", names17[i], sep = "."))))
  }
  normal.theta = do.call(rbind,ll)
  
  #construct initial common trend
  g = function (theta2, theta3, t, kt) {
    d = data.frame (kt, t)
    sm.t = theta3 * t + theta2 # time adjustment
    sm = locpol (kt~t, d, kernel = EpaK, xeval = sm.t) 
    mu1 = sm$lpFit[, 2]
    mu = ts( mu1, start = sm.t[1], frequency = 1)
    return (mu)
  }
  
  for (i in 1: loop2)
  {
    nam15 = paste ("g", i, sep = "")
    assign (nam15, g (normal.theta[i, 2],
                      normal.theta[i, 3], time(kt[,i]),
                      kt[,i]))
  }
  reference.kt = rowMeans (merge, na.rm = TRUE)
  #reference.kt = ts (reference.kt0, start = Sweden$year[1], frequency = 1)
  return ( list (normal.theta, reference.kt))
}