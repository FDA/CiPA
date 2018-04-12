# File:         compute_MCMC_CI.R
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  R script to compute credible intervals from MCMC results.
#               For help, run this script with command line option "-h".
#

#--- specify command line arguments
library(optparse)
parser<-OptionParser()
parser<-add_option(parser, c("-d", "--drug"), default="dofetilide,cisapride,bepridil,verapamil,terfenadine,ranolazine,sotalol,mexiletine,quinidine,ondansetron,diltiazem,chlorpromazine", help="Drug name(s), comma separated [default is 12 CiPA training drugs]")
parser<-add_option(parser, c("-a", "--alpha"), default=0.05, help="Significance level (0-1)")

args<-parse_args(parser)

#--- process arguments
drugstr<-gsub(" ","",args$drug)
drugnames<-strsplit(drugstr, ",")[[1]]
alpha<-args$alpha

for(drug in drugnames){
    outdir<-sprintf("results/%s/",drug)
    print(outdir)

    # get optimal params and MCMC samples
    odf<-read.csv(paste0(outdir,"IC50_optimal.csv"))
    sampdf<-read.csv(paste0(outdir,"IC50_samples.csv"))

    # log transform
    islog<-sapply(colnames(odf), function(col) grepl("IC50$",col))
    colnames(odf)[islog]<-sapply(colnames(odf)[islog], function(param) sub("IC50","log10IC50",param))
    colnames(sampdf)<-colnames(odf)
    notna<-islog & !is.na(odf[1,])
    odf[,notna]<-log10(odf[,notna])
    sampdf[,notna]<-log10(sampdf[,notna])

    cidf<-data.frame()
    for(param in colnames(odf)){
        thisval<-odf[,param]

        # add results using quantile
        if(is.na(thisval))
            out<-c(lower=NA, upper=NA)
        else{
            out<-quantile(sampdf[,param], probs=c(0.025,0.975), na.rm=TRUE)
            names(out)<-c("lower","upper")
        }
        out<-data.frame(param=param, value=thisval, t(out))
        cidf<-rbind(cidf,out)
    }

    cidf$channel<-sapply(cidf$param, function(label) sub("_h$","",sub("_log10IC50$","",label)))
    cidf$param<-mapply(function(label,channel) sub(paste0(channel,"_"),"",label), cidf$param, cidf$channel)

    write.csv(cidf, sprintf("%s/MCMC_CI%g.csv",outdir,100*(1-alpha)), row.names=F, quote=F)
} # for drug

print(warnings())
