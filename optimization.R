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
    out = optim (theta, loss, gr = NULL, t, kt, t.reference, kt.reference, control = list (maxit = 1000))
    loss.mse = out$value
    optimal.theta = out$par
    error.theta = mean ((theta - optimal.theta)^2)
return ( list( loss.mse, optimal.theta, error.theta ))
}