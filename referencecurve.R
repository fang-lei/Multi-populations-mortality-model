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

referencecurve = function (theta2, theta3, kt, kt.null) {
    t = time (kt)  
    sm.t = theta3 * t + theta2 # time adjustment
    t.reference = Sweden$year
    ## common grid for kt and shifted kt.reference
    tmin=max(min(t.reference), min(sm.t))
    tmax=min(max(t.reference), max(sm.t)) 
    i0=which(t.reference>=tmin & t.reference<=tmax)
    if (length(i0)>0) {
    t0=t.reference[i0] 
    d = data.frame (kt, sm.t)
    
    # way 1 - locpol: locpol
    # sm = locpol (kt~sm.t, d, kernel = EpaK, xeval = t0, bw = bw.default) 
    # mu1 = sm$lpFit[, 2]
    
    # way 2 - locpol with Cross Validation bandwidth
    #temp117 = regCVBwSelC(sm.t, kt, deg = 1, kernel = EpaK, interval = c(0, 10))
    #sm = locpol (kt~sm.t, d, kernel = EpaK, xeval = t0, bw = temp117) # time-adjusted kt based on smoothed reference curve
    # mu1 = sm$lpFit[, 2]
    
    # way 3 - sm: sm.regression with optimal smoothing parameter
    h.optimal6 = h.select(sm.t, kt)
    sm = sm.regression(sm.t, kt, h = h.optimal6, eval.points = t0, model = "none", poly.index = 1, display="none")
    mu1 = sm$estimate
    
    # way 4 - KernSmooth: locpoly
    #sm = locpoly (sm.t, kt, kernel = EpaK, bandwidth = bw.default, gridsize = length(t0), range.x = range(t0))
    #mu1 = sm$y
    
    # way 5 - stats: ksmooth
#     sm = ksmooth (sm.t, kt, kernel = "normal", bandwidth = bw.default, n.points = length(t0), range.x = range(t0))
#     mu1 = sm$y
    
    
    mu = ts( mu1, start = t0[1], frequency = 1)}
    else {mu = kt.null}
    return (mu)
}