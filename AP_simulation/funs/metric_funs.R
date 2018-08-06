# File:         metric_funs.R
# Author:       Kelly Chang
# Date:         Sep 2017
# Version:      1.0
# 
# Description:  Helper R functions to find various points in a trace.
#

find_dxdtMax<-function(t, x){
    dxdt<-diff(x) / diff(t)
    i_max<-which.max(dxdt)[1]
    t_max<-t[i_max]
    up_max<-dxdt[i_max]
    list(value=up_max, t=t_max)
}

find_max<-function(t, x){
    i_max<-which.max(x)[1]
    t_max<-t[i_max]
    xmax<-x[i_max]
    list(value=xmax, t=t_max)
}

find_rep<-function(t, x, xrest, idx_max, pct){
    xrep<-xrest+(x[idx_max]-xrest)*(1-pct)
    idx_rep<-which(x<=xrep & t>t[idx_max])[1]
    t_rep<-t[idx_rep]
    list(value=xrep, t=t_rep)
}
