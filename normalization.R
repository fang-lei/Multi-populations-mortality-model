# ------------------------------------------------------------------------------
# Project:     Mortality Model for Multip-populations: A Semiparametric 
#              Comparison Approach
# ------------------------------------------------------------------------------
# Quantlet:    normalization.R
# ------------------------------------------------------------------------------
# Description: Normalize optimal theta.
# ------------------------------------------------------------------------------
# Keywords:    normalization, mortality, Lee-Carter method, commom trend
# ------------------------------------------------------------------------------
# See also:    twopop.R, multipop.R, data.R, optimization.R, referencecurve.R
# ------------------------------------------------------------------------------
# Author:      Lei Fang
# ------------------------------------------------------------------------------

normalization = function (theta) {
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
  return ( normal.theta )
}