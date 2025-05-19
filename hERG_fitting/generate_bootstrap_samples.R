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
library(boot)
library(parallel)

parser<-OptionParser()
parser<-add_option(parser, c("-d", "--drug"), default="quinidine", help="Drug name")
parser<-add_option(parser, c("-f", "--file"), default=NULL, type="character", help="Path to CSV file [required]")
parser<-add_option(parser, c("-n", "--nboots"), default=2000, type="integer", help="Number of bootstrap samples")
parser<-add_option(parser, c("-s", "--seed"), default=100, type="integer", help="Random seed [default 100]")

args<-parse_args(parser)
print(sessionInfo())

# Check required arguments
if(is.null(args$drug)) stop("Missing drug argument!")
if(is.null(args$file)) stop("Missing file path argument!")

#--- process arguments
drug<-args$drug
datafile<-args$file
nboots<-args$nboots
seednum<-args$seed

# Verify file exists
if(!file.exists(datafile)) stop(sprintf("File %s does not exist!", datafile))

#--- sample with replacement using boot package
set.seed(seednum, kind="L'Ecuyer-CMRG")
s<-.Random.seed

    # read in Milnes cell data
    datadf<-read.csv(datafile)
    datadf<-datadf[with(datadf,order(conc,exp,sweep,time)),]

    # data frame of experiments
    expdf<-unique(datadf[,c("conc","exp")])

    # generate bootstrap samples
    boot.out<-boot(data=expdf, statistic=function(datadf, idx){mean(datadf[idx,"exp"])}, R=nboots, strata=as.factor(expdf$conc))

    # save boot object
    outdir<-sprintf("results/%s/",drug)
    system(paste0("mkdir -p ",outdir))
    saveRDS(boot.out,paste0(outdir,"boot_out.rds"))
