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
  mu.ts = ts(mu, start = sm.t[1], frequency = 1)
  return (list(mse, mu.ts))
}
