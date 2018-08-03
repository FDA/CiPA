# File:         hERG_fitting.R
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  R script to fit hERG drug binding model parameters to Milnes
#               protocol data. The drug data to be fitted must be specified
#               with "-d DRUG", where "DRUG" indicates the path to the
#               fractional block data:
#                   data/DRUG.csv
#               For help with other options, run this script with command
#               line option "-h".
#

proc_start<-proc.time()

options(warn=1)

library(optparse)

isWindows<-Sys.info()[["sysname"]]=="Windows"

#--- specify command line arguments
parser<-OptionParser()
parser<-add_option(parser, c("-d", "--drug"), type="character", help="Drug name [required]")
parser<-add_option(parser, c("-i", "--sample"), default="0", help="Bootstrap sample number, range, and/or list of samples to be fitted [default 0 fits to original data]")
parser<-add_option(parser, c("-s", "--seed"), default=100, type="integer", help="Random seed [default 100]")
parser<-add_option(parser, c("-c", "--cores"), default=1, type="integer", help="Number of cores to use during fitting [default 1]")
parser<-add_option(parser, c("-f", "--forking"), default=FALSE, action="store_true", help="Flag to turn on forking for parallelization (not supported in Windows)")
parser<-add_option(parser, c("-l", "--lambda"), type="integer", help="Population size to use for CMA-ES [default 4+floor(3*log(N))]")
parser<-add_option(parser, c("-m", "--maxiter"), type="integer", help="Maximum number of generations for CMA-ES [default 100*N^2]")
parser<-add_option(parser, c("-t", "--tol"), type="double", help="Stopping tolerance for CMA-ES [default 0.5e-12]")

args<-parse_args(parser)

#--- load libraries
library(cmaes)
library(parallel)
library(deSolve)
print(sessionInfo())

#--- required argument
if(is.null(args$drug)) stop("Missing drug argument!")
drug<-args$drug

#--- parse bootstrap list, ranges
bootstr<-args$sample
boot_split<-strsplit(bootstr, ",")[[1]]
boot_list<-lapply(boot_split,
                  function(s){
                      tmp<-as.numeric(strsplit(s, "-")[[1]])
                      if(length(tmp)==0 || length(tmp)>2 || any(is.na(tmp)))
                          stop("Bad sample specifications!")
                      boot_start<-tmp[1]
                      boot_stop<-ifelse(length(tmp)==2, tmp[2], tmp[1])
                      if(boot_start>boot_stop)
                          stop("Invalid sample range!")
                      seq(boot_start, boot_stop) # step size 1
                  })
boot_range<-do.call(c,boot_list)

#--- optional arguments
seednum<-args$seed
cores<-args$cores
usefork<-args$forking
POP_SIZE<-args$lambda
MAX_GENERATION<-args$maxiter
STOPTOL<-args$tol

if(usefork && isWindows){
    print("Windows system detected! Forking not supported, using sockets instead...")
    usefork<-FALSE
}

# cmaes hyperparameters
ctl_list<-list(vectorized=TRUE)
if(!is.null(POP_SIZE)){
    print(sprintf("Using population size: %g",POP_SIZE))
    ctl_list[["lambda"]]<-POP_SIZE
}
if(!is.null(MAX_GENERATION)){
    print(sprintf("Using number of generations: %g",MAX_GENERATION))
    ctl_list[["maxit"]]<-MAX_GENERATION
}
if(!is.null(STOPTOL)){
    print(sprintf("Using stopping tolerance: %g",STOPTOL))
    ctl_list[["stop.tolx"]]<-STOPTOL
}

#--- setup serial or parallel evaluation
setupfile<-"setup_hERG_fitting.R"
usesocket<-FALSE
if(cores>1){
    if(usefork){
        lapplyfun<-function(X, FUN) mclapply(X, FUN, mc.cores=cores, mc.preschedule=FALSE)
    }else{
        cl<-makeCluster(cores)
        invisible(clusterEvalQ(cl,library(deSolve)))
        clusterExport(cl,c("setupfile","drug","isWindows"))

        lapplyfun<-function(X, FUN) clusterApply(cl, X, FUN)
        usesocket<-TRUE
    }
}else{
    lapplyfun<-lapply
}

#--- run range of bootstraps
for(boot_num in boot_range){
    print(sprintf("Running %s bootstrap sample %d on %d cores", drug, boot_num, cores))

    #--- set seed for reproducibility
    set.seed(seednum, kind="L'Ecuyer-CMRG")

    #---setup script
    source(setupfile)
    source("funs/deWrapper.R")

    #--- initial parameter values
    if(boot_num>0){
        pardf<-read.table(sprintf("results/%s/pars.txt",drug)) # start from optimal values
        initval<-setNames(pardf[,2],pardf[,1])
        initpar<-encode_pars(initval[pnames]) # pnames from setupfile
    }else{
        initpar<-runif(length(pnames))*pmax # pnames,pmax from setupfile
    }

    #--- output directory
    outdir<-paste0("results/",drug,"/")
    if(boot_num>0)
        outdir<-paste0(outdir,sprintf("boot/%05d/",boot_num))
    system(paste0("mkdir -p ",outdir))
    print(outdir)

    #--- update setup on cluster
    if(usesocket){
        clusterExport(cl,"boot_num")
        invisible(clusterEvalQ(cl, source(setupfile)))
    }

    #--- prepare vectorized objective function
    genfile<-paste0(outdir,"CMAESgeneration")
    errfile<-paste0(outdir,"CMAESerror")
    if(file.exists(genfile)) file.remove(genfile);
    if(file.exists(errfile)) file.remove(errfile);
    objfun_vec<-function(pop){
        nchr<-ncol(pop)
        chrpernode<-ceiling(nchr/cores)
        nnodes<-min(cores,nchr)
        print(Sys.time())
        errorlist<-lapplyfun(1:nnodes, deWrapper(chrpernode,pop,objfun)) # objfun defined in setupfile
        errormat<-do.call(rbind, errorlist)
        imin<-which.min(errormat[,2])
        write.table(data.frame(t(decode_pars(pop[,imin]))), genfile, sep=" ", row.names=F, col.names=F, append=T)
        write(errormat[imin,2], errfile, ncolumns=1, sep=" ", append=T)
        errormat[,2]
    }

    #--- run CMA-ES
    res<-cma_es(initpar,objfun_vec,lower=0,upper=pmax,control=ctl_list) # pmax from setupfile
    str(res)

    #--- save best results
    pars<-signif(decode_pars(res$par), digits=4)
    names(pars)<-pnames
    write.table(pars, file=paste0(outdir,"pars.txt"), row.names=T, col.names=F, quote=F)

}# for boot_num

if(usesocket) stopCluster(cl)

print(proc.time()-proc_start)
