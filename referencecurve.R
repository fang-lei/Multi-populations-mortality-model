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
    sm.t = theta3 * t + theta2 # time adjustment
    t.reference = Sweden$year
    ## common grid for kt and shifted kt.reference
    tmin=max(min(t.reference), min(sm.t))
    tmax=min(max(t.reference), max(sm.t)) 
    i0=which(t.reference>=tmin & t.reference<=tmax)
    t0=t.reference[i0] 
    d = data.frame (kt, sm.t)
    sm = locpol (kt~sm.t, d, kernel = EpaK, xeval = t0) 
    mu1 = sm$lpFit[, 2]
    mu = ts( mu1, start = t0[1], frequency = 1)
    return (mu)
}