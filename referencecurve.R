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

optimization = function (theta, kt) {
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
  
  return ( list (normal.theta, reference.kt))
}