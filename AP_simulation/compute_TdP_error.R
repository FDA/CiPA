# File:         compute_TdP_error.R
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  R script to perform Torsade de Pointes (TdP) risk
#               classification using ordinal logistic regression and
#               leave-one-out cross validation (LOOCV).
#               For help, run this script with command line option "-h".
#

#--- specify command line arguments
library(optparse)
parser<-OptionParser()
parser<-add_option(parser, c("-t", "--tdpfile"), default="data/CiPA_training_drugs.csv", help="Filepath to table of known TdP risk levels [default data/CiPA_training_drugs.csv]")
parser<-add_option(parser, c("-u", "--uncertainty"), default=FALSE, action="store_true", help="Flag to use simulations with uncertainty inputs for training and cross-validation")
parser<-add_option(parser, c("-o", "--omit"), default=FALSE, action="store_true", help="Flag to use drug simulations with ionic current effects omitted")

args<-parse_args(parser)

#--- load libraries
library(rms)
library(ggplot2)
print(sessionInfo())

#--- get arguments
tdpfile<-args$tdpfile
useUQ<-args$uncertainty
omit<-args$omit

# setup output directory
input<-ifelse(useUQ, "uncertainty", "fixed")
simdir<-sprintf("results/%s/",input)
if(omit){
    outdir_vec<-Sys.glob(paste0(simdir,"drop_*/"))
}else{
    outdir_vec<-simdir
}

#--- function to try ordinal logistic regression
try_lrm<-function(datadf, tol=1e-10, maxit=1e6){
    try({ lrm(CiPA~qNet, data=datadf, penalty=0, x=TRUE, y=TRUE, tol=tol, maxit=maxit) })
}

#--- read in drug TdP risk
drugtable<-read.csv(tdpfile)
drugnames<-as.character(drugtable$drug)

#--- do logistic regression
for(outdir in outdir_vec){
    # read in dataset
    infile<-paste0(outdir,"metrics.rds")
    df<-readRDS(infile)
    df<-df[df$drug!="control" & df$dose!=0,]
    df<-df[!is.na(df$max_dv),] # depolarization failures
    df<-merge(df, drugtable[,c("drug","CiPA")], by="drug", all.x=T)
    df$drug<-factor(df$drug,levels=drugnames)
    df<-df[order(df$drug,df$sample),]
    df$CiPA<-ordered(df$CiPA)

    # perform logistic regression and leave-one-out cross validation for each dose
    errdf<-data.frame()
    probdf<-data.frame()
    cverrdf<-data.frame()
    cvprobdf<-data.frame()
    for(dose in sort(unique(df$dose))){
        for(drug in c(NA,drugnames)){

            # fit logistic regression model
            if(is.na(drug)){
                print(sprintf("training on dose = %g Cmax",dose))
                datadf<-df[df$dose==dose, c("drug","CiPA","qNet")]
                testdf<-datadf
                traindf<-datadf
            }else{
                print(sprintf("cross validating with %s",drug))
                testdf<-datadf[datadf$drug==drug,]
                traindf<-datadf[datadf$drug!=drug,]
                if(nrow(testdf)==0){
                    print("no data to test! skipping...")
                    next
                }
            }

            lmod<-try_lrm(traindf)
            print(lmod)
            if(inherits(lmod, "try-error")){
                if(is.na(drug)){
                    print(sprintf("fitting dose %g failed! skipping...",dose))
                    break
                }else{
                    print(sprintf("cross validating with %s failed! skipping...",drug))
                    next
                }
            }

            # save coefficients
            print(sprintf("Convergence failure: %s",lmod$fail))
            cf<-coefficients(lmod)
            cfvec<-c()
            for(kint in 1:(length(cf)-1))
                cfvec[[paste0("intercept",kint)]]<-cf[[kint]]
            cfvec[["slope"]]<-cf[[length(cf)]]

            # get training/prediction error
            if(is.na(drug)){
                probs<-predict(lmod, type="fitted.ind")
            }else{
                probs<-matrix(predict(lmod, newdata=testdf, type="fitted.ind"), ncol=length(levels(y0)))
            }
            y0<-testdf$CiPA
            yPred<-apply(probs,1, function(x) which.max(x))
            pred_err<-mean(abs(yPred-as.integer(y0)))

            # append to data frame
            if(is.na(drug)){
                newrow<-data.frame(dose=dose, t(cfvec), error=pred_err)
                errdf<-rbind(errdf,newrow)
            }else{
                newrow<-data.frame(drug=drug, dose=dose, t(cfvec), error=pred_err)
                cverrdf<-rbind(cverrdf,newrow)
            }

            # detailed results
            newcols<-sapply(levels(y0), function(s) paste0("predict_",s), USE.NAMES=F)
            predictions<-factor(newcols[yPred], levels=newcols)
            tmpdf<-as.data.frame(tapply(1:nrow(testdf), list(drug=as.character(testdf$drug), predict=predictions), FUN=length))
            tmpdf[is.na(tmpdf)]<-0
            tmpdf$drug<-factor(rownames(tmpdf), levels=drugnames)
            tmpdf$dose<-dose
            tmpdf<-merge(drugtable[,c("drug","CiPA")],tmpdf)
            tmpdf<-tmpdf[order(tmpdf$drug),]
            rownames(tmpdf)<-NULL
            if(is.na(drug)){
                probdf<-rbind(probdf,tmpdf)
            }else{
                cvprobdf<-rbind(cvprobdf,tmpdf)
            }

        } # for drug
    } # for dose

    outfile<-paste0(outdir,"training_probs.csv")
    write.csv(probdf, outfile, row.names=F, quote=F)
    print(head(probdf))

    outfile<-paste0(outdir,"training_errors.csv")
    write.csv(errdf, outfile, row.names=F, quote=F)
    print(head(errdf))

    outfile<-paste0(outdir,"LOOCV_probs.csv")
    write.csv(cvprobdf, outfile, row.names=F, quote=F)
    print(head(cvprobdf))

    outfile<-paste0(outdir,"LOOCV_errors.csv")
    write.csv(cverrdf, outfile, row.names=F, quote=F)
    print(head(cverrdf))
}# for outdir

#--- get training results to plot
df<-data.frame()
for(outdir in outdir_vec){
    infile<-paste0(outdir,"training_errors.csv")
    tmpdf<-read.csv(infile)

    if(omit){
        dropcurrent<-sub("/","",sub(".*drop_","",outdir))
        tmpdf$drop<-dropcurrent
    }else{
        tmpdf$drop<-"none"
    }

    df<-rbind(df,tmpdf)
}# for outdir

system(paste0("mkdir -p figs"))
dropstr<-ifelse(omit, "dropcurrent_", "")
figfile<-sprintf("figs/%s%s_training_errors.pdf",dropstr,input)
print(figfile)
pdf(figfile, width=8, height=4.5)
p<-ggplot(df, aes(x=dose, y=error, color=drop))
p<-p+geom_point(size=2)+geom_line(size=0.75)
p<-p+ylab("Training error")
p<-p+xlab(~"Concentration (\u00D7"~C[max]*")")
p<-p+theme_bw()
if(!omit)
    p<-p+theme(legend.position="none")
print(p)
dev.off()

#--- get LOOCV results to plot
df<-data.frame()
for(outdir in outdir_vec){
    infile<-paste0(outdir,"LOOCV_errors.csv")
    tmpdf<-read.csv(infile)

    if(omit){
        dropcurrent<-sub("/","",sub(".*drop_","",outdir))
        tmpdf$drop<-dropcurrent
    }else{
        tmpdf$drop<-"none"
    }

    df<-rbind(df,tmpdf)
}# for outdir

system(paste0("mkdir -p figs"))
dropstr<-ifelse(omit, "dropcurrent_", "")
figfile<-sprintf("figs/%s%s_LOOCV_errors.pdf",dropstr,input)
print(figfile)
pdf(figfile, width=8, height=4.5)
p<-ggplot(df, aes(x=dose, y=error, color=drop))
p<-p+stat_summary(fatten=3, fun.data=mean_sdl)
p<-p+ylab("Prediction error")
p<-p+xlab(~"Concentration (\u00D7"~C[max]*")")
p<-p+theme_bw()
if(!omit)
    p<-p+theme(legend.position="none")
print(p)
dev.off()
