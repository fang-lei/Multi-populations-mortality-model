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
# See also:    twopop.R, multipop.R, data.R, referencecurve.R, normalization.R
# ------------------------------------------------------------------------------
# Author:      Lei Fang, Juhyun Park
# ------------------------------------------------------------------------------

optimization = function (theta, kt, kt.reference) {
  t = time (kt)
  t.reference = time (kt.reference)
  ### loss function
    loss = function (theta, t, kt, t.reference, kt.reference) {
      ## assume theta[1]>0, theta[3]>0, if not, take absolute values
      theta1 = abs(theta[1])
      theta2 = theta[2]
      theta3 = abs(theta[3])
     
      if (theta[1]<0 | theta[3]<0) warning("theta1 or theta3 <0, abs value is used")
      sm.t = (t.reference - theta2) /theta3 # time adjustment
      ## common domain for kt and time-adjusted kt.reference
      tmin=max(min(t), min(sm.t))
      tmax=min(max(t), max(sm.t)) 
      i0=which(t>=tmin & t<=tmax) # index for common domain
      if (length(i0)>0){
        t0=t[i0]  
        ## smooth interpolation of shifted kt.reference on common grid    
        dref=data.frame(sm.t=sm.t, kt.reference=kt.reference)
        
        # way 1 - locpol: locpol
        #sm = locpol(kt.reference~sm.t, dref, kernel = EpaK, xeval=t0, bw = bw.default) # time-adjusted kt based on smoothed reference curve
        #mu = theta1 * sm$lpFit[,2] # modelled kt based on theta
        
        # way 2 - locpol with Cross Validation bandwidth
        #temp117 = regCVBwSelC(sm.t,kt.reference, deg = 1, kernel = EpaK, interval = c(0, 10))
        #sm = locpol(kt.reference~sm.t, dref, kernel = EpaK, xeval=t0, bw = temp117) # time-adjusted kt based on smoothed reference curve
        #mu = theta1 * sm$lpFit[,2] # modelled kt based on theta
        
        # way 3 - sm: sm.regression with optimal smoothing parameter
        h.optimal3 = h.select(sm.t, kt.reference)
        sm = sm.regression(sm.t, kt.reference, h = h.optimal3, eval.points = t0, model = "none", poly.index = 1, display="none")
        mu = theta1 * sm$estimate
        
        # way 4 - KernSmooth: locpoly
        #sm = locpoly (sm.t, kt.reference, kernel = EpaK, bandwidth = bw.default, gridsize = length(t0), range.x = range(t0))
        #mu = theta1 * sm$y
        
        # way 5 - stats: ksmooth
        # sm = ksmooth (sm.t, kt.reference, kernel = "normal", bandwidth = bw.default, n.points = length(t0), range.x = range(t0))
        # mu = theta1 * sm$y
        
        ## mean squared error at common grid points
        mse = mean ((kt[i0]-mu)^2) # mse of kt and the modelled one
      } else { mse=1e9 }
      return (mse)
    }
  ### parameter estimation
    conv=1
    theta0 = theta
    while(conv!=0) { # check convergence
      out = optim (theta0, loss, gr = NULL, t, kt, t.reference, kt.reference, control = list (maxit = 1000)) # optimization
      conv=out$convergence
      theta0=out$par
    }  
    result=c(out$par, out$value, out$convergence)
    
  ### check results on graph
    theta=result[1:4]
    theta1=abs(theta[1])
    theta2=theta[2]
    theta3=abs(theta[3])
    
    sm.t = (t.reference - theta2) /theta3 # time adjustment
    ## common grid for kt and shifted kt.reference
    tmin=max(min(t), min(sm.t))
    tmax=min(max(t), max(sm.t)) 
    i0=which(t>=tmin & t<=tmax)
    t0=t[i0]  
    ## shifted curves (kt.hat)
    dref=data.frame(sm.t=sm.t, kt.reference=kt.reference)
    
    # way 1 - locpol: locpol
    #sm = locpol(kt.reference~sm.t, dref, kernel = EpaK, xeval=sm.t, bw = bw.default) # time-adjusted kt based on smoothed reference curve
    #mu = theta1 * sm$lpFit[,2]  # modelled kt based on theta
    
    # way 2 - locpol with Cross Validation bandwidth
    #temp117 = regCVBwSelC(sm.t,kt.reference, deg = 1, kernel = EpaK, interval = c(0, 10))
    #sm = locpol(kt.reference~sm.t, dref, kernel = EpaK, xeval=sm.t, bw = temp117) # time-adjusted kt based on smoothed reference curve
    #mu = theta1 * sm$lpFit[,2]  # modelled kt based on theta
    
    # way 3 - sm: sm.regression with optimal smoothing parameter
    h.optimal4 = h.select (sm.t, kt.reference)
    sm = sm.regression(sm.t, kt.reference, h = h.optimal4, eval.points = sm.t, model = "none", poly.index = 1, display="none")
    mu = theta1 * sm$estimate
    
    # way 4 - KernSmooth: locpoly
    #sm = locpoly (sm.t, kt.reference, kernel = EpaK, bandwidth = bw.default, gridsize = length(sm.t), range.x = range(sm.t))
    #mu = theta1 * sm$y
    
    # way 5 - stats: ksmooth
#     sm = ksmooth (sm.t, kt.reference, kernel = "normal", bandwidth = bw.default, n.points = length(sm.t), range.x = range(sm.t))
#     mu = theta1 * sm$y
    
    
    mu = ts (mu, start = sm.t[1], frequency = theta3)
    
    
    # way 1 - locpol: locpol
    #sm0 = locpol(kt.reference~sm.t, dref, kernel = EpaK, xeval=t0, bw = bw.default) # time-adjusted kt based on smoothed reference curve
    #mu0 = theta1 * sm0$lpFit[,2]  # modelled kt based on theta
    
    # way 2 - locpol with Cross Validation bandwidth
    #temp117 = regCVBwSelC(sm.t,kt.reference, deg = 1, kernel = EpaK, interval = c(0, 10))
    #sm0 = locpol(kt.reference~sm.t, dref, kernel = EpaK, xeval=t0, bw = temp117) # time-adjusted kt based on smoothed reference curve   
    #mu0 = theta1 * sm0$lpFit[,2]  # modelled kt based on theta
    
    # way 3 - sm: sm.regression with optimal smoothing parameter
    h.optimal5 = h.select (sm.t, kt.reference)
    sm0 = sm.regression(sm.t, kt.reference, h = h.optimal5, eval.points = t0, model = "none", poly.index = 1, display="none")
    mu0 = theta1 * sm0$estimate
    
    # way 4 - KernSmooth: locpoly
    #sm0 = locpoly (sm.t, kt.reference, kernel = EpaK, bandwidth = bw.default, gridsize = length(t0), range.x = range(t0))
    #mu0 = theta1 * sm0$y
    
    # way 5 - stats: ksmooth
#     sm0 = ksmooth (sm.t, kt.reference, kernel = "normal", bandwidth = bw.default, n.points = length(t0), range.x = range(t0))
#     mu0 = theta1 * sm0$y
    
    
    mu0 = ts (mu0, start = t0[1], frequency = 1)
    
return ( list( result, mu, mu0,tmin,tmax ))
}