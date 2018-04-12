# File:         combine_results.R
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  R script to aggregate action potential (AP) simulation
#               results.
#               For help, run this script with command line option "-h".
#

#--- specify command line arguments
library(optparse)
parser<-OptionParser()
parser<-add_option(parser, c("-n", "--num_samples"), default=0, type="integer", help="How many uncertainty samples to combine [default 0 combines fixed-input results]")
parser<-add_option(parser, c("-o", "--omit"), default=FALSE, action="store_true", help="Flag to combine drug simulations with ionic current effects omitted")

args<-parse_args(parser)

#--- get arguments
nsamp<-args$num_samples
omit<-args$omit

if(nsamp<0) stop("Number of samples must be nonnegative!")

# setup output directory
input<-ifelse(nsamp==0, "fixed", "uncertainty")
simdir<-sprintf("results/%s/",input)
if(omit){
    outdir_vec<-Sys.glob(paste0(simdir,"drop_*/"))
}else{
    outdir_vec<-simdir
}

# get simulation results
for(outdir in outdir_vec){
    filename<-"metrics.rds"
    if(input=="fixed"){
        df_list<-list()
        for(infile in Sys.glob(paste0(outdir,"*/",filename))){
            if(!file.exists(infile)) next
            drug<-sub("/.*","",sub(outdir,"",infile))
            df<-readRDS(infile)
            df$drug<-drug
            df$sample<-NA
            df_list[[drug]]<-df
        }
        df<-do.call(rbind,df_list)
        rownames(df)<-NULL
    }else{
        df_list<-list()
        for(drugdir in Sys.glob(paste0(outdir,"*/"))){
            drug<-sub("/","",sub(outdir,"",drugdir))
            samp_list<-list()
            for(isamp in 1:nsamp){
                infile<-sprintf("%s/%05d/%s",drugdir,isamp,filename)
                if(!file.exists(infile)) next
                df<-readRDS(infile)
                df$drug<-drug
                df$sample<-isamp
                samp_list[[isamp]]<-df
            }
            df_list[[drug]]<-do.call(rbind,samp_list)
            rownames(df_list[[drug]])<-NULL
        }
        df<-do.call(rbind,df_list)
        rownames(df)<-NULL
    }
    print(head(df))

    outfile<-paste0(outdir,filename)
    saveRDS(df, outfile)
}# for outdir
