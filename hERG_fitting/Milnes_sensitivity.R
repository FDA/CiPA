# File:         Milnes_sensitivity.R
# Author:       Kelly Chang
# Date:         Sep 2017
# Version:      1.0
# 
# Description:  R script to do sensitivity analysis using hERG bootstrap
#               results for the 12 CiPA training drugs. 
#               New drugs can be specified by running this script with the
#               command line option "-d DRUG".
#               For help with other options, run this script with command
#               line option "-h".
#

proc_start<-proc.time()

#--- specify command line arguments
library(optparse)

parser<-OptionParser()
parser<-add_option(parser, c("-d", "--drug"), default="dofetilide,cisapride,bepridil,verapamil,terfenadine,ranolazine,sotalol,mexiletine,quinidine,ondansetron,diltiazem,chlorpromazine", help="Drug name(s), comma separated [default is 12 CiPA training drugs]")

args<-parse_args(parser)

#--- load libraries
library(deSolve)
library(FME)
library(ggplot2)
print(sessionInfo())

isWindows<-Sys.info()[["sysname"]]=="Windows"

#--- process arguments
drugstr<-gsub(" ","",args$drug)
drugnames<-strsplit(drugstr, ",")[[1]]

#--- function for running Milnes simulation
solveMilnes<-function(fitpars, conc, times=peaktimes){
    # run simulation
    drugsweeps<-hergmodel$run_drug(fitpars, conc) # from setupfile
    if(length(drugsweeps)==0)
        stop("Solving error!")

    outdf<-data.frame()
    startt<-0
    for(i in 1:nbeats){
        ctltime<-ctlsweeps[[i]][,"time"]
        idxout<-ctltime%in%times
        drugtime<-ctltime[idxout]+startt
        drugO<-drugsweeps[[i]][idxout,"O"]
        drugFrac<-drugO/ctlsweeps[[i]][idxout,"O"]
        outdf<-rbind(outdf,data.frame(time=drugtime,O=drugO,frac=drugFrac))
        startt<-startt+tail(fulltimes,1)
    } # for beats

    outdf
} # solveMilnes

for(drug in drugnames){
    #--- setupfile used in fitting
    print(drug)
    boot_num<-0
    setupfile<-"setup_hERG_fitting.R"
    source(setupfile)
    conc_vec<-as.numeric(concvec) # concvec from setupfile is characters
    print(conc_vec)

    #--- get fitting and bootstrap results
    outdir<-sprintf("results/%s/",drug)
    tmp<-read.table(paste0(outdir,"pars.txt"), header=FALSE, row.names=1)
    fitpars<-t(tmp)[1,]
    bootpars<-na.omit(read.csv(paste0(outdir,"boot_pars.csv")))

    bootpars<-bootpars[,pnames] # pnames from setupfile
    fitpars<-fitpars[pnames]

    print(fitpars)
    print(head(bootpars))

    #--- get sensitivity of model output to variability in parameters
    datadf<-data.frame()
    sensdf<-data.frame()
    for(conc in conc_vec){
        # output times from data (1st sweep)
        deptime<-fracdata[[as.character(conc)]][[1]][,"time"]

        # specify drug conc
        print(sprintf("concentration = %g nM",conc))
        sR<-sensRange(solveMilnes, parms=fitpars, sensvar="frac", parInput=bootpars, num=nrow(bootpars), conc=conc, times=deptime)

        # save summary of traces
        summ.sR<-summary(sR)

        # additional quantiles
        summ.sR[,"q02.5"]<-apply(sR[,(length(fitpars)+1):ncol(sR)], 2, FUN=function(x) quantile(x, probs = 0.025))
        summ.sR[,"q97.5"]<-apply(sR[,(length(fitpars)+1):ncol(sR)], 2, FUN=function(x) quantile(x, probs = 0.975))

        # only save depolarization time
        for(i in 1:nbeats){
            beatstart<-(i-1)*tail(fulltimes,1)
            beatlim<-range(deptime)+beatstart
            preddf<-as.data.frame(summ.sR[summ.sR[,"x"]>=beatlim[1] & summ.sR[,"x"]<=beatlim[2],])
            rownames(preddf)<-NULL
            preddf$x<-preddf$x-beatstart
            preddf$conc<-conc
            preddf$sweep<-i
            sensdf<-rbind(sensdf,preddf)

            # get data
            meandf<-fracdata[[as.character(conc)]][[i]]
            meandf<-meandf[meandf$time%in%preddf$x,]
            ilow<-unique(c(seq(1,nrow(meandf),by=5), nrow(meandf))) # every 5th point
            meandf<-meandf[ilow,]
            meandf$conc<-conc
            meandf$sweep<-i
            datadf<-rbind(datadf,meandf)
        }
    }

    # save results
    saveRDS(sensdf, paste0(outdir,"sensRange_summary.rds"))

    # plot results
    sensdf$x<-sensdf$x-940
    datadf$time<-datadf$time-940

    ustr<-" nM"
    if(min(sensdf$conc)>10){
        ustr<-" \u03BCM"
        sensdf$conc<-sensdf$conc/1e3
        datadf$conc<-datadf$conc/1e3
    }
    clevels<-sort(unique(sensdf$conc), decreasing=TRUE)
    clabels<-sapply(clevels, function(conc) paste0(conc,ustr))
    sensdf$conc<-factor(sensdf$conc, levels=clevels)
    datadf$conc<-factor(datadf$conc, levels=clevels)

    p<-ggplot(sensdf, aes(x=x/1000, group=conc))
    p<-p+facet_wrap(~sweep, labeller="label_both")
    p<-p+geom_ribbon(aes(ymin=q02.5, ymax=q97.5, fill=conc), alpha=0.3)
    p<-p+geom_line(data=datadf, aes(x=time/1000, y=frac, color=conc), size=0.25, alpha=0.5)
    p<-p+coord_cartesian(ylim=c(-0.2,1.2))
    p<-p+scale_fill_manual(NULL, values=1:length(clevels), labels=clabels)
    p<-p+scale_color_manual(NULL, values=1:length(clevels), labels=clabels)
    p<-p+scale_x_continuous(breaks=c(0,5,10))
    p<-p+scale_y_continuous(breaks=c(0,0.5,1))
    p<-p+xlab("Time (s)")
    p<-p+ylab("Fractional current")
    p<-p+guides(fill=guide_legend(ncol=2), color=guide_legend(ncol=2))
    p<-p+theme_bw()
    p<-p+theme(plot.title=element_text(hjust=0.5), legend.position=c(0.875,0.15), 
            panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    ggsave(paste0("figs/",drug,"_sensRange.png"), p, width=8, height=4)
}# for drug

print(proc.time()-proc_start)
