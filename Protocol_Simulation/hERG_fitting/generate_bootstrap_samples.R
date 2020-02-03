# File:         generate_bootstrap_samples.R
# Author:       Kelly Chang
# Date:         Sep 2017
# Version:      1.0
# 
# Description:  R script to generate bootstrap samples from hERG channel
#               fractional block data for 12 CiPA training drugs.
#               New drugs can be specified by running this script with the
#               command line option "-d DRUG".
#               The required input file(s) are comma-separated value (CSV)
#               files located at:
#                   data/DRUG.csv
#               For help with other options, run this script with command
#               line option "-h".
#

#--- specify command line arguments
library(optparse)

myDrugDefault = "astemizole,azimilide,clarithromycin,clozapine,disopyramide,domperidone,droperidol,ibutilide,loratadine,metoprolol,nifedipine,nitrendipine,pimozide,risperidone,tamoxifen,vandetanib"

parser<-OptionParser()
parser<-add_option(parser, c("-d", "--drug"), default=myDrugDefault, help="Drug name(s), comma separated [default is 12 CiPA training drugs]")
parser<-add_option(parser, c("-n", "--nboots"), default=2000, type="integer", help="Number of bootstrap samples")
parser<-add_option(parser, c("-s", "--seed"), default=100, type="integer", help="Random seed [default 100]")

args<-parse_args(parser)

#--- load libraries
library(boot)
library(parallel)
print(sessionInfo())

#--- process arguments
drugstr<-gsub(" ","",args$drug)
drugnames<-strsplit(drugstr, ",")[[1]]
nboots<-args$nboots
seednum<-args$seed

#--- sample with replacement using boot package
set.seed(seednum, kind="L'Ecuyer-CMRG")
s<-.Random.seed
for(tid in seq_along(drugnames)){
    drug<-drugnames[tid]
    print(drug)

    # set random seed
    if(tid>1) s<-nextRNGStream(s)
    .Random.seed<-s

    # read in Milnes cell data
    datadf<-read.csv(paste0("data/",drug,".csv"))
    datadf<-datadf[with(datadf,order(conc,exp,sweep,time)),]

    # data frame of experiments
    expdf<-unique(datadf[,c("conc","exp")])

    # generate bootstrap samples
    boot.out<-boot(data=expdf, statistic=function(datadf, idx){mean(datadf[idx,"exp"])}, R=nboots, strata=as.factor(expdf$conc))

    # save boot object
    outdir<-sprintf("results/%s/",drug)
    system(paste0("mkdir -p ",outdir))
    saveRDS(boot.out,paste0(outdir,"boot_out.rds"))
}# for tid
