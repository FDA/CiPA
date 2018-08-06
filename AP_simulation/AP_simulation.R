# File:         AP_simulation.R
# Author:       Kelly Chang
# Date:         Sep 2017
# Version:      1.0
# 
# Description:  R script to run action potential (AP) simulations with the
#               IKr-dynamic ORd model. Drug effects are simulated using
#               previously fitted hERG binding kinetic parameters and Hill
#               equation parameters for six other ionic currents. The drug
#               parameters can be either fixed-input (i.e. optimal), or
#               uncertainty-input parameters selected from their sampling
#               distributions.  Several AP and calcium transient metrics are
#               saved, as well as the qNet metric described by Dutta et al.:
#               https://www.frontiersin.org/article/10.3389/fphys.2017.00616
#               For help with simulation options, run this script with
#               command line option "-h".
#

proc_start<-proc.time()

options(warn=1)

#--- specify command line arguments
library(optparse)
parser<-OptionParser()
parser<-add_option(parser, c("-d", "--drug"), default="control", help="Drug name")
parser<-add_option(parser, c("-i", "--sample"), default="0", help="Sample number, range, and/or list [default 0 for fixed parameters]")
parser<-add_option(parser, c("-x", "--dose"), default="0", help="Dose as a multiple of therapeutic concentration [default 0]")
parser<-add_option(parser, c("-o", "--omit"), default=FALSE, action="store_true", help="Flag to simulate drug with effects on hERG, ICaL, INaL, INa, IKs, Ito, and IK1 omitted, one at a time")
parser<-add_option(parser, c("-n", "--no_opt_scaling"), default=FALSE, action="store_true", help="Flag to turn off optimized ionic conductance scaling factors")
parser<-add_option(parser, c("-r", "--dropcurrent"), help="Simulate drug with effects on one ionic current omitted [hERG, ICaL, INaL, INa, IKs, Ito, or IK1]")
parser<-add_option(parser, c("-c", "--cmaxfile"), default="data/CiPA_training_drugs.csv", help="Filepath to table of therapeutic concentrations [default data/CiPA_training_drugs.csv]")
parser<-add_option(parser, c("-b", "--hergpath"), default="../hERG_fitting/results/", help="Path to hERG fitting results [default ../hERG_fitting/results/]")
parser<-add_option(parser, c("-m", "--hillpath"), default="../Hill_fitting/results/", help="Path to Hill fitting results [default ../Hill_fitting/results/]")

args<-parse_args(parser)

#--- load libraries and helper functions
library(deSolve)
source("funs/get_pars.R")
source("funs/metric_funs.R")

print(sessionInfo())

isWindows<-Sys.info()[["sysname"]]=="Windows"

#--- get arguments
drug<-tolower(args$drug)
sampstr<-args$sample
dosestr<-args$dose
omit<-args$omit
use_opt_scaling<-!args$no_opt_scaling
dropstr<-args$dropcurrent
cmaxfile<-args$cmaxfile
bootpath<-args$hergpath
mcmcpath<-args$hillpath

#--- parse sample list, ranges
samp_split<-strsplit(sampstr, ",")[[1]]
samp_list<-lapply(samp_split,
                  function(s){
                      tmp<-as.numeric(strsplit(s, "-")[[1]])
                      if(length(tmp)==0 || length(tmp)>2 || any(is.na(tmp)))
                          stop("Bad sample specifications!")
                      samp_start<-tmp[1]
                      samp_stop<-ifelse(length(tmp)==2, tmp[2], tmp[1])
                      if(samp_start>samp_stop)
                          stop("Invalid sample range!")
                      seq(samp_start, samp_stop) # step size 1
                  })
samples<-do.call(c,samp_list)
if(any(samples<0)) stop("Sample number must be nonnegative!")

#--- parse dose list, ranges
if(drug=="control"){
    doses<-0
}else{
    dose_split<-strsplit(dosestr, ",")[[1]]
    dose_list<-lapply(dose_split,
                      function(s){
                          tmp<-as.numeric(strsplit(s, "-")[[1]])
                          if(length(tmp)==0 || length(tmp)>2 || any(is.na(tmp)))
                              stop("Bad dose specifications!")
                          dose_start<-tmp[1]
                          dose_stop<-ifelse(length(tmp)==2, tmp[2], tmp[1])
                          if(dose_start>dose_stop)
                              stop("Invalid dose range!")
                          seq(dose_start, dose_stop) # step size 1
                      })
    doses<-do.call(c,dose_list)
}
if(any(doses<0)) stop("Dose must be nonnegative!")

#--- parse dropcurrent
all_currents<-c("hERG","ICaL","INaL","INa","IKs","Ito","IK1")
drop_vec<-NA
if(omit){
    drop_vec<-all_currents
}else if(!is.null(dropstr)){
    if(!dropstr%in%all_currents)
        stop(sprintf("Unrecognized current: %s!",dropstr))
    drop_vec<-dropstr
}

#--- load model
extension<-ifelse(isWindows, ".dll", ".so")
dyn.load(paste0("models/newordherg_qNet",extension))

#--- get therapeutic concentration (Cmax)
drugtable<-read.csv(cmaxfile)
cmax<-drugtable[tolower(as.character(drugtable$drug))==drug,"therapeutic"] # should be in nanomolar
if(length(cmax)==0){
    if(drug=="control"){
        cmax<-0
    }else{
        cmax<-1
        print(sprintf("Cmax undefined, interpreting dose as nanomolar concentrations."))
    }
}else if(length(cmax)==1){
    print(sprintf("Cmax set to %g nM, interpreting dose as multiples of Cmax.",cmax))
}else{
    stop("Multiple entries for %s therapeutic concentration!",drug)
}

#--- run simulations, optionally removing drug effects on one current at a time
CL<-2000
beats<-1000
nsave<-250
print(sprintf("CL = %s ms, analyzing %d out of %d beats", CL, nsave, beats))

for(isamp in samples){
    print(sprintf("Running %s, sample %d, doses %s", drug, isamp,
                  ifelse(length(doses)==1,doses,args$dose)))

    for(dropcurrent in drop_vec){

        #--- setup output directory
        if(drug=="control" || isamp==0){
            outdir<-"results/fixed/"
        }else{
            outdir<-"results/uncertainty/"
        }
        if(!is.na(dropcurrent)){
            print(sprintf("Omitting drug effects on %s...", dropcurrent))
            outdir<-paste0(outdir,"drop_",dropcurrent,"/")
        }
        outdir<-paste0(outdir,drug,"/")
        if(isamp>0) outdir<-paste0(outdir,sprintf("%05d/",isamp))
        system(paste0("mkdir -p ",outdir))

        #--- setup simulation outputs
        outfile<-paste0(outdir,"metrics.rds")
        outdf<-data.frame()

        statefile<-paste0(outdir,"states_max_dv.rds")
        statesdf<-data.frame()

        chkpt_intv<-250
        checkfile<-paste0(outdir,"checkpoints.rds")
        checksdf<-data.frame()

        #--- run simulations
        for(dose in doses){
            conc<-dose*cmax
            print(sprintf("dose = %sx Cmax, concentration = %g nM...", dose, conc))

            pars<-get_pars("models/newordherg_pars.txt", drug, conc, isamp, bootpath, mcmcpath, dropcurrent, use_opt_scaling)
            print(pars)

            tmp<-read.table(sprintf("models/newordherg_states_CL%d.txt",CL), col.names=c("param","value"))
            states<-setNames(tmp$value, tmp$param)
            states[["D"]]<-conc

            # add qNet variable to initial states
            states<-c(states, qNet=0)

            fulltimes<-0:CL

            beatstates<-list()
            beatdf<-data.frame()
            params<-c("INa","INaL","Ito","ICaL","IKr","IKs","IK1","dv")
            for(i in 1:beats){
                # run simulation
                try({out <- ode(states, fulltimes, "derivs", pars,dllname="newordherg_qNet",initfunc="initmod", nout=length(params),rtol=1e-6,atol=1e-6,method="lsoda")}); 
                if(!exists("out")||inherits(out,"try-error")||length(out[,1])!=length(fulltimes)
                   || !all(out[,1]==fulltimes) || any(is.nan(out))){
                    stop("solving error!")
                }
                colnames(out)[1+length(states)+seq_along(params)]<-params

                # save checkpoint
                if(i>1 && i%%chkpt_intv==1)
                    checksdf<-rbind(checksdf,data.frame(drug=drug, CL=CL, dose=dose, beat=i, t(states)))

                # save metrics
                if(i>beats-nsave){
                    beatstates[[i-(beats-nsave)]]<-data.frame(beat=i, t(states)) # save all beats temporarily
                    tmpdf<-data.frame(beat=i, t(out[nrow(out),"qNet"]))

                    # default values
                    catd50<-NA
                    catd90<-NA
                    apd50<-NA
                    apd90<-NA

                    # Mirams metric calculations
                    t<-out[,"time"]
                    v<-out[,"v"]
                    ca<-out[,"cai"]
                    t_start<-t[length(t)] # start looking for EAD
                    t_end<-t[length(t)] # stop looking for EAD

                    # transmembrane potential
                    pt_dvdtmax<-find_dxdtMax(t, v)
                    t_up<-pt_dvdtmax[["t"]]
                    tmpdf$dvdtmax<-pt_dvdtmax[["value"]]

                    pt_max<-find_max(t, v)
                    t_max<-pt_max[["t"]]
                    vmax<-pt_max[["value"]]
                    i_max<-which(t==t_max)
                    tmpdf$vmax<-vmax
                    if(vmax>=0) t_start<-t_max

                    vrest<-min(v[t<t_max])
                    tmpdf$vrest<-vrest
                    if(vmax>=0 && vrest< -40){ # repolarization
                        t_rep<-find_rep(t, v, vrest, i_max, 0.3)[["t"]]
                        if(!is.na(t_rep)){
                            t_start<-t_rep # start looking for EAD

                            t_rep<-find_rep(t, v, vrest, i_max, 0.5)[["t"]]
                            if(!is.na(t_rep)){
                                apd50<-t_rep-t_up

                                t_rep<-find_rep(t, v, vrest, i_max, 0.9)[["t"]]
                                if(!is.na(t_rep)){
                                    apd90<-t_rep-t_up
                                    t_end<-t_rep # stop looking for EAD
                                }
                            }
                        }
                    }
                    tmpdf$APD50<-apd50
                    tmpdf$APD90<-apd90

                    tmpdf$max_dv<-ifelse(t_end-t_start>1,
                                         max(out[t>t_start & t<t_end,"dv"]),
                                         NA)

                    # intracellular calcium
                    pt_max<-find_max(t, ca)
                    t_max<-pt_max[["t"]]
                    camax<-pt_max[["value"]]
                    i_max<-which(t==t_max)
                    tmpdf$camax<-camax

                    carest<-min(ca[t<t_max])
                    tmpdf$carest<-carest

                    caup<-carest+(camax-carest)*0.9
                    t_up<-t[which(ca>caup)[1]]

                    t_rep<-find_rep(t, ca, carest, i_max, 0.5)[["t"]]
                    if(!is.na(t_rep)){
                        catd50<-t_rep-t_up

                        t_rep<-find_rep(t, ca, carest, i_max, 0.9)[["t"]]
                        if(!is.na(t_rep)){
                            catd90<-t_rep-t_up
                            t_end<-t_rep
                        }
                    }
                    tmpdf$CaTD50<-catd50
                    tmpdf$CaTD90<-catd90

                    # append metrics
                    beatdf<-rbind(beatdf,tmpdf)
                }

                # set state variables for next beat
                pidx<-match(names(states),colnames(out),nomatch=0)
                states[pidx!=0]<-out[dim(out)[1],pidx]
                states[length(states)]<-0 # reset integration

            } # for beats

            # save info from last beat, or beat with max_dv
            max_vec<-beatdf$max_dv
            namaxdv<-is.na(max_vec)
            EADs<-sum(max_vec[!namaxdv]>0)
            noDepols<-sum(namaxdv)
            noRepols<-sum(is.na(beatdf$APD90))
            if(noDepols==length(max_vec)){
                keep<-nrow(beatdf)
                print(sprintf("depolarization failure, keeping beat %d",beatdf[keep,"beat"]))
            }else{
                keep<-which.max(max_vec)
                print(sprintf("max_dv on beat %d",beatdf[keep,"beat"]))
            }

            saveRDS(checksdf,checkfile)

            statesdf<-rbind(statesdf,data.frame(dose=dose, beatstates[[keep]]))
            saveRDS(statesdf,statefile)

            outdf<-rbind(outdf,data.frame(dose=dose, beatdf[keep,], EADs=EADs, noDepols=noDepols, noRepols=noRepols))
            saveRDS(outdf,outfile)

        } # for dose
    } # for dropcurrent
} # for sample

print(proc.time()-proc_start)
