# File:         IC50_mcmc.R
# Author:       Kelly C. Chang
#                Zhihua Li
#               based on code by Jose Vicente
# Date:         Oct 2017
# Version:      1.1
# 
# Description:  R script to perform Markov-chain Monte Carlo (MCMC)
#               simulation for Hill equations modeling drug block of
#               ionic currents. The drug to be fitted must be specified with
#               "-d DRUG". The data for the specified drug should be located
#               in "data/drug_block.csv" (default), or at the file path
#               specified by "-f FILEPATH".
#               For help with other options, run this script with command
#               line option "-h".
#

proc_start<-proc.time()

options(warn=1)

library(optparse)

#--- specify command line arguments
parser<-OptionParser()
parser<-add_option(parser, c("-d", "--drug"), type="character", help="Drug name [required]")
parser<-add_option(parser, c("-f", "--filepath"), default="data/drug_block.csv", help="Path to file containing patch clamp fractional block data for the specified drug [default data/drug_block.csv].")
parser<-add_option(parser, c("-s", "--seed"), default=100, type="integer", help="Random seed [default 100]")
parser<-add_option(parser, c("-n", "--nsamples"), default=2000, type="integer", help="Number of samples to collect [default 2000]")
parser<-add_option(parser, c("-b", "--burnin"), default=10000, type="integer", help="Burnin length [default 10000]")
parser<-add_option(parser, c("-t", "--thin"), default=10, type="integer", help="Thinning rate [default 10]")
parser<-add_option(parser, c("-c", "--channel"), type="character", default="all", help="Channel or Current name")
args<-parse_args(parser)

#--- load libraries
library(FME)
library(coda)
print(sessionInfo())

#--- required argument
if(is.null(args$drug)) stop("Missing drug argument!")
drug<-args$drug

#--- optional arguments
datafile<-args$filepath
seednum<-args$seed
nsamp<-args$nsamples
burnin<-args$burnin
thin<-args$thin
desiredchannels<-args$channel

sigdigits=4
addburn<-burnin # increase burnin by addburn if convergence test fails

# parameter bounds for log(IC50) and Hill coefficient (note IC50 in nM)
lowBds<-c(log(1e-10), 0)
uppBds<-c(log(1e10), 10)

# Hill equation residuals
HillRes<-function(params, datadf){
    datadf$block - Hillfun(params, datadf$conc)$block
}

# Hill equation
Hillfun<-function(params, conc){
    data.frame(conc=conc, block=100*(1-1/(1+(conc/exp(params[1]))^params[2])))
}

# Hill equation self start function
# note: logIC50 is used to keep IC50 from becoming negative
SSHill<-selfStart(~ 100*(1-1/(1+(conc/exp(logIC50))^h)),
                  function(mCall, data, LHS){
                      xy <- sortedXyData(mCall[["conc"]], LHS, data)
                      if(nrow(xy) < 4)
                          stop("Too few distinct x values to fit Hill equation")
                      pars<-c(logIC50=NA, h=NA)
                      eps<-0.01
                      if(all(xy[["y"]]<eps)){ # 0% block for all conc
                          pars[["logIC50"]]<-log(max(xy[["x"]]))
                          pars[["h"]]<-1
                      }else if(all(xy[["y"]]>100-eps)){ # 100% block for all conc, should handle with care
                          pars[["logIC50"]]<-log(min(xy[["x"]]))
                          pars[["h"]]<-1
                      }else{ # linear fit
                          ignore<- xy[["y"]]<eps | xy[["y"]]>100-eps
                          xx<-log(xy[["x"]][!ignore])
                          yy<-1/(1-xy[["y"]][!ignore]/100)-1
                          yy<-log(yy)
                          xxyy<-data.frame(xx=xx, yy=yy)
                          lm.out<-lm(yy ~ xx, data=xxyy)
                          cf<-coef(lm.out)
                          pars[["logIC50"]]<- -cf[[1]]/cf[[2]]
                          pars[["h"]]<-cf[[2]]
                      }
                      return(pars)
                  }, c("logIC50","h"))

# read in patch clamp data
datadf<-read.csv(datafile)
datadf<-datadf[datadf$drug==drug,]
datadf<-datadf[with(datadf,order(drug,conc,channel)),]

channels<-sort(unique(datadf$channel))

# save plots
figdir<-"figs/"
system(paste0("mkdir -p ",figdir))
pdf(paste0(figdir,drug,"_nls_mcmc.pdf"),width=6,height=4,useDingbats=FALSE)

# save samples to table
drugdir<-sprintf("results/%s/",drug)
system(paste0("mkdir -p ",drugdir))

# fit Hill equation for each channel
opt_list<-list()
col_list<-list()
if(desiredchannels == "all")
     desiredchannels <- channels
for(channel in desiredchannels){
    tstr<-sprintf("%s, %s channel",drug,channel)
    print(tstr)

    # get data frame of cells for this experiment
    expdf<-datadf[datadf$drug==drug & datadf$channel==channel,]
    expdf<-expdf[!is.na(expdf$block),] # remove NA

    # prepare variables to store output
    col_list[[paste0(channel,"_IC50")]]<-rep(NA,nsamp)
    col_list[[paste0(channel,"_h")]]<-rep(NA,nsamp)
    opt_list[[paste0(channel,"_IC50")]]<-NA
    opt_list[[paste0(channel,"_h")]]<-NA
   
    # fit model using modFit (nls.lm)
    initpar<-getInitial(block ~ SSHill(conc, logIC50, h), expdf)
    meanconc<-mean(expdf$conc)
    ic50seeds<-c(meanconc,meanconc*0.5,meanconc*0.3,meanconc*0.75,0.00001, 0.0001, 0.001, 0.01, 0.1, 1,100,300,500,1000,3000,5000,10000,30000,50000,100000, 300000, 500000,1000000,3000000,5000000)
    trystarts<-cbind(log(ic50seeds),rep(0.9,length(ic50seeds)))
    trystarts<-rbind(initpar,trystarts)
    tryi<-1
    while(tryi==1 || (inherits(tryout, "try-error") && tryi<=nrow(trystarts)) ){
        tryout<-try({
            #optim is more robust; always use it first to narrow it down unless ...
           if(!any(is.na(trystarts[tryi,]))){
            realinitpars <- optim(trystarts[tryi,],function(p) sum(HillRes(p,expdf)^2),lower=lowBds, upper=uppBds,method="L-BFGS-B")$par
            }else{                              #if init is NA, still need to go through modFit to get tryout
            realinitpars<-c(logIC50=1,h=1)                #but change init to a random pair
            }
            mf<-modFit(f=HillRes, p=realinitpars, datadf=expdf, 
                       lower=lowBds, upper=uppBds, method="Marq")

            smf<-summary(mf)
            print(smf)
        })
        tryi<-tryi+1
    }

    # plot data only or do MCMC
    xvals<-range(expdf$conc)
    if(inherits(tryout, "try-error")){
        print("Fitting error! Skipping MCMC...")

        # plot data
        par(mfrow=c(1,1))
        plot(expdf$conc, expdf$block, log="x",  
             main=tstr, xlab="Concentration (nM)", ylab="Block (%)", 
             xlim=range(xvals), ylim=c(0,100))
        next
    }

    # variables used to plot dose-response curve
    plot_IC50<-exp(mf$par[1])
    plot_h<-mf$par[2]
    xvals<-c(xvals, plot_IC50/exp(2), plot_IC50*exp(2))
    plot_x<-exp(seq(from=min(log(xvals)), to=max(log(xvals)), length.out=100))
    plot_y<-SSHill(plot_x, log(plot_IC50), plot_h)

    # save IC50 values
    opt_list[[paste0(channel,"_IC50")]]<-signif(plot_IC50, digits=sigdigits)
    opt_list[[paste0(channel,"_h")]]<-signif(plot_h, digits=sigdigits)

    # initialize with modFit results, as recommended in FME documentation
    startp<-setNames(mf$par, c("logIC50","h"))
    Covar<-smf$cov.scaled * 2.4^2/2
    s2prior<-smf$modVariance

    # set seed for reproducibility
    set.seed(seednum)
    burnin<-addburn
    converged<-FALSE
    while(!converged){
        # run MCMC
        tryout<-try({
            MCMC<-modMCMC(f=HillRes, p=startp, datadf=expdf, 
                          lower=lowBds, upper=uppBds,
                          jump=Covar, var0=s2prior, wvar0=1,
                          burninlength=burnin, niter=burnin+nsamp*thin, outputlength=nsamp,
                          updatecov=10, ntrydr=2)
            print(summary(MCMC))
        })
        if(inherits(tryout, "try-error") || MCMC$naccepted<nsamp){ #MCMC failed; let's try changing jump 
            print("Error running MCMC! Changing jump to default value")

            tryout<-try({
            MCMC<-modMCMC(f=HillRes, p=startp, datadf=expdf, 
                          lower=lowBds, upper=uppBds,
                          jump=NULL, var0=s2prior, wvar0=1,
                          burninlength=burnin, niter=burnin+nsamp*thin, outputlength=nsamp,
                          updatecov=10, ntrydr=2)
            print(summary(MCMC))
        })
        }
          if(inherits(tryout, "try-error") || MCMC$naccepted<nsamp){ #MCMC run still failed
            print("Error running MCMC! Removing all IC50 values...")

            # Remove previously fitted optimal values
            opt_list[[paste0(channel,"_IC50")]]<-NA
            opt_list[[paste0(channel,"_h")]]<-NA
            break
        }

        # check for convergence
        gd<-geweke.diag(as.mcmc(MCMC$pars))
        if(all(abs(gd[["z"]])<qnorm(0.975))){
            converged<-TRUE
            saveRDS(MCMC, paste0(drugdir,sprintf("%s_mcmc.rds",channel)))

            # get sensitivity
            sR<-sensRange(func=Hillfun, parms=mf$par, parInput=MCMC$pars, conc=plot_x, num=nsamp)
            saveRDS(sR, paste0(drugdir,sprintf("%s_sensRange.rds",channel)))

            # plot sensitivity
            par(mfrow=c(1,1))
            plot(summary(sR), quant=TRUE, obs=expdf[,c("conc","block")],
                 log="x", xlab="Concentration (nM)", ylab="Block (%)", ylim=c(0,100), main=tstr)
            lines(plot_x, plot_y, col="red") # plot modFit

            # plot MCMC results
            par(mfrow=c(1,1))
            plot(MCMC$pars, main=tstr)
            plot(MCMC)
            hist(MCMC)

            # save results to table
            col_list[[paste0(channel,"_IC50")]]<-signif(exp(MCMC$pars[,1]), digits=sigdigits) # transform back to IC50
            col_list[[paste0(channel,"_h")]]<-signif(MCMC$pars[,2], digits=sigdigits)
        }else{
            print("Geweke diagnostic indicates lack of convergence, increasing burnin...")
            set.seed(seednum)
            burnin<-burnin+addburn
        }
    }# while not converged
}# for channel

optdf<-do.call(cbind,opt_list)
print(sprintf("Saving optimal IC50s for %s",drug))
print(head(optdf))
write.csv(optdf, paste0(drugdir,"IC50_optimal.csv"), row.names=F, quote=F)

drugdf<-do.call(cbind,col_list)
print(sprintf("Saving samples for %s",drug))
print(head(drugdf))
write.csv(drugdf, paste0(drugdir,"IC50_samples.csv"), row.names=F, quote=F)

dev.off()

print(proc.time()-proc_start)