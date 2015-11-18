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
      mu = theta1 * sm$lpFit[,2] + theta4 # modelled kt based on initial theta
      mse = mean ((kt-mu)^2) # mse of kt and the modelled one
      return (mse)
    }
    out = optim (theta, loss, gr = NULL, t, kt, t.reference, kt.reference, control = list (maxit = 1000))
    optimal.theta = out$par
    dref = data.frame (kt, t)
    sm.t1 = (t - optimal.theta[2]) /optimal.theta[2] # time adjustment
    sm1 = locpol (kt~t, dref, kernel = EpaK, xeval = sm.t1) # time-adjusted kt based on smoothed reference curve
    shift.kt1 = optimal.theta[1] * sm1$lpFit[,2] + optimal.theta[4] # modelled kt based on initial theta
    shift.kt = ts (shift.kt1, start = sm.t1[1], frequency = 1)
    error.theta = mean ((theta - optimal.theta)^2)
    error.kt = mean ((kt - shift.kt)^2)
return ( list( optimal.theta, shift.kt, error.theta, error.kt ))
}