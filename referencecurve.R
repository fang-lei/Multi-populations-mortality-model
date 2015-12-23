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
# See also:    twopop.R, multipop.R, data.R, optimization.R, normalization.R
# ------------------------------------------------------------------------------
# Author:      Lei Fang
# ------------------------------------------------------------------------------

referencecurve = function (theta2, theta3, kt) {
    t = time (kt)  
    d = data.frame (kt, t)
    sm.t = theta3 * t + theta2 # time adjustment
    ## common grid for kt and shifted kt.reference
    tmin=max(min(t), min(sm.t))
    tmax=min(max(t), max(sm.t)) 
    i0=which(t>=tmin & t<=tmax)
    t0=t[i0]  
    sm = locpol (kt~t, d, kernel = EpaK, xeval = t0) 
    mu1 = sm$lpFit[, 2]
    mu = ts( mu1, start = t0[1], frequency = 1/theta3)
    return (mu)
}