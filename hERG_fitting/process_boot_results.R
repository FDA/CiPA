# File:         process_boot_results.R
# Author:       Kelly Chang
# Date:         Sep 2017
# Version:      1.0
# 
# Description:  R script to do postprocessing of hERG fitting results for
#               the 12 CiPA training drugs. This script aggregates
#               bootstrapping results, plots the bootstrap distributions of
#               the fitted parameters, and computes confidence intervals.
#               New drugs can be specified by running this script with the
#               command line option "-d DRUG".
#               For help with other options, run this script with command
#               line option "-h".
#

#--- specify command line arguments
library(optparse)

parser<-OptionParser()
parser<-add_option(parser, c("-d", "--drug"), default="dofetilide,cisapride,bepridil,verapamil,terfenadine,ranolazine,sotalol,mexiletine,quinidine,ondansetron,diltiazem,chlorpromazine", help="Drug name(s), comma separated [default is 12 CiPA training drugs]")
parser<-add_option(parser, c("-a", "--alpha"), default=0.05, help="Significance level (0-1)")

args<-parse_args(parser)

#--- load libraries
library(boot)
print(sessionInfo())

#--- process arguments
drugstr<-gsub(" ","",args$drug)
drugnames<-strsplit(drugstr, ",")[[1]]
alpha<-args$alpha

sigdigits<-4 # round parameters to significant digits

for(drug in drugnames){
    #--- read in bootstraps
    outdir<-sprintf("results/%s/",drug)
    boot.out<-readRDS(paste0(outdir,"boot_out.rds"))
    nboots<-boot.out$R

    #--- get optimal params (starting point)
    tmp<-read.table(paste0(outdir,"pars.txt"), header=FALSE, row.names=1)
    parStart<-t(tmp)[1,]

    #--- create bootstrap data frame
    df_list<-list()
    for(boot_num in 1:nboots){
        bootdir<-sprintf("results/%s/boot/%05d/", drug, boot_num)
        outfile<-paste0(bootdir,"pars.txt")
        #save.image()
        if(file.exists(outfile)){
         #   print(paste0("file ",outfile," exists!\n"))
            tmp<-read.table(outfile, header=FALSE, col.names=c("param","value"))
          #  print(tmp)
          #  print(paste0("levels(tmp$param) is ",levels(tmp$param)))
          #  levels(tmp$param)[1]<-"Kmax"                #some old files call it Kf still
            
            pars<-setNames(tmp$value, tmp$param)
            newrow<-data.frame(boot=boot_num, t(pars))
        }else{
            newrow<-data.frame(boot=boot_num, t(parStart))
            newrow[,names(parStart)]<-NA
        }
        df_list[[boot_num]]<-newrow
    }
    df<-do.call(rbind, df_list)
    df<-df[order(df$boot),names(parStart)]
    print(tail(df))

    #--- add slope parameter
    parStart[["slope"]]<-signif(parStart[["Kmax"]]/parStart[["halfmax"]], digits=sigdigits)
    df$slope<-signif(df$Kmax/df$halfmax, digits=sigdigits)

    #--- save samples to file
    write.csv(df, paste0(outdir,"boot_pars.csv"), row.names=F, quote=F)

    #--- construct boot object
    boot.out$t0<-parStart
    boot.out$t<-as.matrix(df)

    #--- log transform
    islog<-!names(parStart)%in%c("n","Vhalf")
    names(parStart)[islog]<-sapply(names(parStart)[islog], function(param) paste0("log10",param))
    parStart[islog]<-log10(parStart[islog])
    df[,islog]<-log10(df[,islog])

    #--- plot and calculate quantiles
    figdir<-"figs/"
    system(paste0("mkdir -p ",figdir))
    pdf(paste0(figdir,drug,"_boot.pdf"),width=8,height=4.5)
    cidf<-data.frame()
    for(i in 1:length(parStart)){
        print(parStart[i])
        param<-names(parStart)[i]

        plot(boot.out, index=i)
        title(main=param)

        out<-quantile(df[,i], probs=c(alpha/2,1-alpha/2), na.rm=TRUE)
        names(out)<-c("lower","upper")
        out<-data.frame(param=param, value=parStart[[i]], t(out))
        cidf<-rbind(cidf,out)
    }
    dev.off()
    write.csv(cidf, sprintf("%s/boot_CI%g.csv",outdir,100*(1-alpha)), row.names=F, quote=F)
}# for drug
