# ------------------------------------------------------------------------------
# Project:     Mortality Model for Multip-populations: A Semiparametric 
#              Comparison Approach
# ------------------------------------------------------------------------------
# Quantlet:    optimization.R
# ------------------------------------------------------------------------------
# Description: Optimize theta with kt and reference curve for commom kt.
# ------------------------------------------------------------------------------
# Keywords:    optimization, mortality, Lee-Carter method, smoothing
# ------------------------------------------------------------------------------
# See also:    twopop.R, multipop.R, data.R, referencecurve.R
# ------------------------------------------------------------------------------
# Author:      Lei Fang
# ------------------------------------------------------------------------------

optimization = function (theta, kt, kt.reference) {
  t = time (kt)
  t.reference = time (kt.reference)
    loss = function (theta, t, kt, t.reference, kt.reference) {
      theta1 = theta[1]
      theta2 = theta[2]
      theta3 = theta[3]
      theta4 = theta[4]
      dref = data.frame (kt.reference, t.reference)
      sm.t = (t - theta2) /theta3 # time adjustment
      sm = locpol (kt.reference~t.reference, dref, kernel = EpaK, xeval = sm.t) # time-adjusted kt based on smoothed reference curve
      mu = theta1 * sm$lpFit[,2] + theta4 # modelled kt based on theta
      mse = mean ((kt-mu)^2) # mse of kt and the modelled one
      return (mse)
    }
    out = optim (theta, loss, gr = NULL, t, kt, t.reference, kt.reference, control = list (maxit = 1000)) # optimization
    loss.mse = out$value # loss value from optimization
    optimal.theta = out$par # optimal theta
    error.theta = mean ((theta - optimal.theta)^2) # error (difference / MSE ) between initail theta and optimal theta, for convergence monitor
    # generate the shifted kt based on optiaml theta for plot purpose (comparison of reference curve vs shifted kt)
    dref = data.frame (kt.reference, t.reference)
    sm.t1 = (t - optimal.theta[2]) / optimal.theta[3] # time adjustment
    sm1 = locpol (kt.reference~t.reference, dref, kernel = EpaK, xeval = sm.t1) # time-adjusted kt based on smoothed reference curve
    shift.kt1 = optimal.theta[1] * sm1$lpFit[,2] + optimal.theta[4] # modelled kt based on theta
    shift.kt = ts (shift.kt1, start = sm.t1[1], frequency = 1)
return ( list( loss.mse, optimal.theta, error.theta, shift.kt ))
}