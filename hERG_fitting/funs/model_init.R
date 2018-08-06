# File:         model_init.R
# Author:       Kelly Chang
#               Zhihua Li
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Helper R function to initialize a drug binding model and run
#               voltage clamp simulations at different drug concentrations.
#

model_init<-function(modelname, states, pars, pnames, fulltimes, eventdata, nbeats){
    force(modelname)
    force(nbeats)
    pidx<-match(names(pars), pnames, nomatch=0)
    closingtimes<-seq(0,3*60*1000,10000)
    sweeptimes<-c()
    sweepevents<-NULL

    mymodel<-list()
    mymodel$run_simulation<-function(initstates, pars, timepoints, events=NULL){
        pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]
        try({out <- ode(initstates, timepoints, "derivs", pars, dllname=modelname,
            initfunc="initmod", nout=0, rtol=1e-3, atol=1e-6, method="lsoda",
            events=events)});
        if(!exists("out")||inherits(out,"try-error")||length(out[,1])!=length(timepoints) || !all(out[,1]==timepoints) || any(is.nan(out)))
            return(c())
        out
    }

    mymodel$run_sweeps<-function(initstates, pars){
        states<-initstates
        sweeps<-list()
        for(i in 1:nbeats){
            out<-mymodel$run_simulation(states, pars, sweeptimes, sweepevents)
            if(length(out)==0)
                return(out)
            states[sidx!=0]<-out[nrow(out),sidx]
            sweeps[[i]]<-out
        }
        sweeps
    }

    states["D"]<- 0
    states["V"]<- -85
    out<-mymodel$run_simulation(states, pars, closingtimes)
    if(length(out)==0)
        stop("Solving error during model initialization!")
    sidx<-match(names(states),colnames(out),nomatch=0)
    states[sidx!=0]<-out[nrow(out),sidx]
    ctlstates<-states
    ctlpars<-pars

    sweeptimes<-fulltimes
    sweepevents<-list(data=eventdata)
    ctlsweeps<-mymodel$run_sweeps(ctlstates, ctlpars)
    if(length(ctlsweeps)==0)
        stop("Solving error during sweep initialization!")

    mymodel$controlsweeps<-function() ctlsweeps

    mymodel$run_drug<-function(ind, conc){
        states<-ctlstates
        states["D"]<-conc

        pars<-ctlpars
        pars[pidx!=0]<-ind[pidx]
        pars["timeout"]<-10

        # close channel in presence of drug
        states["V"]<- -80
        pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]
        out<-mymodel$run_simulation(states, pars, closingtimes)
        if(length(out)==0)
            return(out)
        states[sidx!=0]<-out[nrow(out),sidx]

        drugsweeps<<-mymodel$run_sweeps(states, pars)
        drugsweeps
    }

    mymodel
}
