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
      theta4 = theta[4]
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
        sm = locpol(kt.reference~sm.t, dref, kernel = EpaK, xeval=t0) # time-adjusted kt based on smoothed reference curve    
        mu = theta1 * sm$lpFit[,2] + theta4 # modelled kt based on theta
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
    theta4=theta[4]
    sm.t = (t.reference - theta2) /theta3 # time adjustment
    ## common grid for kt and shifted kt.reference
    tmin=max(min(t), min(sm.t))
    tmax=min(max(t), max(sm.t)) 
    i0=which(t>=tmin & t<=tmax)
    t0=t[i0]  
    ## shifted curves (kt.hat)
    dref=data.frame(sm.t=sm.t, kt.reference=kt.reference)
    sm = locpol(kt.reference~sm.t, dref, kernel = EpaK, xeval=sm.t) # time-adjusted kt based on smoothed reference curve    
    mu = theta1 * sm$lpFit[,2] + theta4 # modelled kt based on theta
    mu = ts (mu, start = sm.t[1], frequency = theta3)
    
    sm0 = locpol(kt.reference~sm.t, dref, kernel = EpaK, xeval=t0) # time-adjusted kt based on smoothed reference curve    
    mu0 = theta1 * sm0$lpFit[,2] + theta4 # modelled kt based on theta
    mu0 = ts (mu0, start = t0[1], frequency = 1)
    
return ( list( result, mu, mu0,tmin,tmax ))
}