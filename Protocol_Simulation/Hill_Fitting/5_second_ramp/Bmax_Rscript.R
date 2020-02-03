# File:         IC50_mcmc.R
# Author:       Kelly C. Chang
#                Zhihua Li
#               based on code by Jose Vicente
# Date:         Oct 2017
# Version:      2.0
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
#source("calculate_real_acceptance_rate.R")
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
set.seed(seednum)
# parameter bounds for log(IC50) and Hill coefficient (note IC50 in nM)
lowBds<-c(log(1e-10)*0.99, 0.5*1.01, 70*1.01)
uppBds<-c(log(1e10)*0.99, 2.0*0.99, 100*0.99)

#full bounds include std of the error term
fulllowBds<-c(log(1e-10), 0.5, 70,0)
fulluppBds<-c(log(1e10), 2.0, 100,50)


# Hill equation residuals
HillRes<-function(params, datadf){
	datadf$block - Hillfun(params, datadf$conc)$block
}

# Hill equation
Hillfun<-function(params, conc){
	data.frame(conc=conc, block=params[3]/(1+(exp(params[1])/conc)^params[2]))
}

#likelihood fun
Hilllikelihood<-function(params, datadf){
	#here params has three variables: ic50, h, and std (no.4), and Bmax (no.3)
	res<- HillRes(params[1:3], datadf)
	N<- dim(datadf)[1];  delta<- params[4] 
	loglikelihood<- -N*log(delta)-1/(2*delta^2)*sum(res^2)
	#write.table(data.frame(paste(params,collapse="_"),loglikelihood),"lh",row.names=F,col.names=F,quote=F,append=T)
	-2*(-N*log(delta)-1/(2*delta^2)*sum(res^2))
	
}

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
#myDir <- sprintf("datadf_investigation/%s/%s/",drug,channel)
#dir.create(myDir,recursive = TRUE)
#write.csv(datadf,sprintf("%s/mydata.csv",myDir),row.names = FALSE)
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
	#print(tstr)
	myDir <- sprintf("datadf_investigation/%s/%s/",drug,channel)
	dir.create(myDir,recursive = TRUE)
	write.csv(datadf,sprintf("%s/mydata.csv",myDir),row.names = FALSE)
	
	# get data frame of cells for this experiment
	expdf<-datadf[datadf$drug==drug & datadf$channel==channel,]
	expdf<-expdf[!is.na(expdf$block),] # remove NA
	
	# prepare variables to store output
	col_list[[paste0(channel,"_IC50")]]<-rep(NA,nsamp)
	col_list[[paste0(channel,"_h")]]<-rep(NA,nsamp)
	opt_list[[paste0(channel,"_IC50")]]<-NA
	opt_list[[paste0(channel,"_h")]]<-NA
	
	# fit model using modFit (nls.lm)
	
	meanconc<-mean(expdf$conc)
	ic50seeds<-c(meanconc,meanconc*0.5,meanconc*0.3,meanconc*0.75,0.00001, 0.0001, 0.001, 0.01, 0.1, 1,100,300,500,1000,3000,5000,10000,30000,50000,100000, 300000, 500000,1000000,3000000,5000000)
	trystarts<-cbind(log(ic50seeds),rep(0.9,length(ic50seeds)),rep(0.1*length(ic50seeds)))
	
	colnames(trystarts)<-c("logIC50","h","Bmax")
	tryi<-1
	
	while(tryi==1 || (inherits(tryout, "try-error") && tryi<=nrow(trystarts)) ){
		tryout<-try({
					#optim is more robust; always use it first to narrow it down unless ...
					if(!any(is.na(trystarts[tryi,]))){
						realinitpars <- optim(trystarts[tryi,],function(p) sum(HillRes(p,expdf)^2),lower=lowBds, upper=uppBds,method="L-BFGS-B")
						print(realinitpars)
						realinitpars <- realinitpars$par
					}else{                              #if init is NA, still need to go through modFit to get tryout
						realinitpars<-c(logIC50=1,h=1,Bmax = 100)                #but change init to a random pair
					}
					
					mf<-modFit(f=HillRes, p=realinitpars, datadf=expdf, 
							lower=lowBds, upper=uppBds, method="Marq", hessian = TRUE)
					
				})
		tryi<-tryi+1
	}
	#print("mf hessian:")
	#print(mf$hessian)
	#print("full mf:")
	#print(mf)
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
	plot_h<-mf$par[2];            plot_Bmax<-mf$par[3]
	xvals<-c(xvals, plot_IC50/exp(2), plot_IC50*exp(2))
	#because sometimes around IC50 is not covering all concentrations tested, need to expand more
	xvals<-c(xvals, min(expdf$conc)/10, max(expdf$conc)*10)
	plot_x<-exp(seq(from=min(log(xvals)), to=max(log(xvals)), length.out=100))
	plot_y<-Hillfun(c(log(plot_IC50), plot_h, plot_Bmax),plot_x)$block
	
	# save IC50 values
	opt_list[[paste0(channel,"_IC50")]]<-signif(plot_IC50, digits=sigdigits)
	opt_list[[paste0(channel,"_h")]]<-signif(plot_h, digits=sigdigits)
	opt_list[[paste0(channel,"_Bmax")]]<-signif(plot_Bmax, digits=sigdigits)
	# initialize with modFit results, as recommended in FME documentation
	startp<-setNames(mf$par, c("logIC50","h","Bmax"))
	startp[4]<- sd(mf$residuals)
	
	Covar <- NULL  
	s2prior <- NULL  #this will assume error variance is 1: so MCMC function won't really use error term in calculating likelihood; but I'll do it myself in Hilllikelihood
	weightvar<- 0
	#JUMP <- c(log(1.02),0.02*startp[2],0.02*startp[3])
	
	#Run the MCMC chain for a much shorter number of samples, to try and get a good value for the jump parameter.
#====
	#myJump <- solve(mf$hessian)
	bypass <- FALSE
	
	if(bypass == FALSE)
	{
		
		print("Not bypassing marching strategy.")
	}
	

	
	if(bypass == FALSE)
	{
		print("Calculating jump parameter")
#		burnin.jumpcalc <- 1000
#		nsamp.jumpcalc <- 2000
		jumpFactor <- 1.00
		upperAcceptableLimit <- 0.90
		lowerAcceptableLimit <- 0.30
		acceptanceRatio <- Inf
		failureCounter <- 0
		failureCounterMax <- 200
#		print(sprintf("Burn-In for jump calculation = %d.",burnin.jumpcalc))
#		print(sprintf("# of samples for jump calculation = %d.",nsamp.jumpcalc))
		print(sprintf("Lower acceptable limit = %0.5g.",lowerAcceptableLimit))
		print(sprintf("Upper acceptable limit = %0.5g.",upperAcceptableLimit))
		print(sprintf("Starting; failure counter = %d. Maximum allowed failures = %d.",failureCounter,failureCounterMax))
		
		converged <- FALSE
		burnin <- args$burnin
		converged <- FALSE
		while(((acceptanceRatio < lowerAcceptableLimit | acceptanceRatio > upperAcceptableLimit) & (failureCounter <= failureCounterMax)) | (!converged))
		{
			
			print(sprintf("Inputting jump factor = %0.5g",jumpFactor))
			
			myChain <- modMCMC(f = Hilllikelihood,
					p = startp,
					datadf = expdf, 
					lower = fulllowBds,
					upper = fulluppBds,
					jump = jumpFactor * startp,
					var0 = s2prior,
					wvar0 = weightvar,
					burninlength = burnin,
					niter = burnin + nsamp * thin,
					outputlength = nsamp,
					updatecov=10,
					ntrydr = 2)
			
			nAccepted <- myChain$naccepted
			acceptanceRatio <- nAccepted/(burnin + nsamp * thin)
			#acceptanceRatio <- calculate_real_acceptance_rate(mcmc_chain = myChain, chain_length = nsamp)
			
			print(sprintf("Acceptance ratio = %0.5g",acceptanceRatio))
			if(acceptanceRatio > upperAcceptableLimit)
			{
				failureCounter <- failureCounter + 1
				print(sprintf("Failures = %d.",failureCounter))
				print("acceptanceRatio too high; increasing jump factor.")
				jumpFactor <- jumpFactor * 1.05
			}
			
			if(acceptanceRatio < lowerAcceptableLimit)
			{
				failureCounter <- failureCounter + 1
				print(sprintf("Failures = %d.",failureCounter))
				print("acceptanceRatio too low; decreasing jump factor.")
				jumpFactor <- jumpFactor * 0.95
			}
			
			if((acceptanceRatio < upperAcceptableLimit) & (acceptanceRatio > lowerAcceptableLimit))
			{
				print("Successfully found an estimate for the jump parameter.")
				print("Now testing if the Geweke diagnostic is OK.")
				geweke_diagnostic_test <- geweke.diag(as.mcmc(myChain$pars))
				
				if(all(abs(geweke_diagnostic_test[["z"]]) < qnorm(0.975)))
				{
					print("Geweke diagnostic: OK!")
					converged <- TRUE
				}
				else
				{
					print("Geweke diagnostic: NOT OK! Increasing burnin...")
					burnin <- burnin + addburn
				}
				print(sprintf("Using jumpFactor = %0.5g",jumpFactor))
			}
			
			if(failureCounter >= failureCounterMax)
			{
				stop("Could not find a good initial jump factor.")
			}
			
			
			
			print(sprintf("Jump factor = %0.5g.",jumpFactor))
		}
		
		
	}
	print("myChain:")
	print(myChain)
	print("Starting parameters:")
	print(startp)
	
#save.image()
	# set seed for reproducibility
	
	#burnin<-addburn
	#converged<-FALSE
	
	#Test to make sure same results are obtained as were obtained when the jump-finding code terminated.
#	print("RUNNING TEST SECTION:")
#	print("JUMP FACTOR:")
#	print(jumpFactor)
#	myChain_test <- modMCMC(f = Hilllikelihood,
#			p = startp,
#			datadf = expdf, 
#			lower = fulllowBds,
#			upper = fulluppBds,
#			jump = jumpFactor * startp,
#			var0 = s2prior,
#			wvar0 = weightvar,
#			burninlength = burnin.jumpcalc,
#			niter = burnin.jumpcalc + nsamp.jumpcalc,
#			updatecov=10,
#			ntrydr = 2)
#	print("BEGIN TEST RUN RESULTS:")
#	print(summary(myChain_test))
#	print("END TEST RUN RESULTS.")
	stuff <- FALSE
	while(!stuff){
		stuff <- TRUE
		# run MCMC
#		tryout<-try({
#					MCMC<-modMCMC(f = Hilllikelihood,
#							p = startp,
#							datadf = expdf, 
#							lower = fulllowBds,
#							upper = fulluppBds,
#							jump = jumpFactor*startp,
#							var0 = s2prior,
#							wvar0 = weightvar,
#							burninlength = burnin,
#							niter = burnin + nsamp*thin,
#							outputlength = nsamp,
#							updatecov = 10,
#							ntrydr = 2)
#					MCMC <- myChain
#					
#					print(summary(MCMC))
#				})
#		
#		if(inherits(tryout, "try-error") || MCMC$naccepted<nsamp){ #MCMC run still failed
#			print("Error running MCMC! Removing all IC50 values...")
#			
#			# Remove previously fitted optimal values
#			opt_list[[paste0(channel,"_IC50")]]<-NA
#			opt_list[[paste0(channel,"_h")]]<-NA
#			break
#		}
		MCMC <- myChain
		# check for convergence
		gd<-geweke.diag(as.mcmc(MCMC$pars))
		if(all(abs(gd[["z"]])<qnorm(0.975))){
			converged<-TRUE
			saveRDS(MCMC, paste0(drugdir,sprintf("%s_mcmc.rds",channel)))
			
#save.image()
			# get sensitivity
			sR<-sensRange(func=Hillfun, parms=mf$par, parInput=MCMC$pars[,1:2], conc=plot_x, num=nsamp)
			saveRDS(sR, paste0(drugdir,sprintf("%s_sensRange.rds",channel)))
			
			# plot sensitivity
			par(mfrow=c(1,1))
			plot(summary(sR), quant=TRUE, obs=expdf[,c("conc","block")],
					log="x", xlab="Concentration (nM)", ylab="Block (%)", ylim=c(-50,150), main=tstr)
			lines(plot_x, plot_y, col="red") # plot modFit
			
			# plot MCMC results
			par(mfrow=c(1,1))
			plot(MCMC$pars, main=tstr)
			plot(MCMC)
			hist(MCMC)
			
			# save results to table
	col_list[[paste0(channel,"_IC50")]]<-signif(exp(MCMC$pars[,1]), digits=sigdigits) # transform back to IC50
	col_list[[paste0(channel,"_h")]]<-signif(MCMC$pars[,2], digits=sigdigits)
	col_list[[paste0(channel,"_Bmax")]]<-signif(MCMC$pars[,3], digits=sigdigits)
		}else{
			print("Geweke diagnostic indicates lack of convergence, increasing burnin...")
			set.seed(seednum)
			burnin <- burnin + addburn
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
#save.image()
print(proc.time()-proc_start)