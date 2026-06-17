# File:         compute_TdP_error.R
# Author:       Kelly Chang
#                Zhihua Li
# Date:         Nov 2017
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
parser<-add_option(parser, c("-t", "--tdpfile"), default="../nanion8rand_AP_simulation/data/newCiPA.csv", help="Filepath to drug list")
parser<-add_option(parser, c("-u", "--uncertainty"), default=TRUE, action="store_true", help="Flag to use simulations with uncertainty inputs for training and cross-validation")
parser<-add_option(parser, c("-o", "--omit"), default=FALSE, action="store_true", help="Flag to use drug simulations with ionic current effects omitted")
parser<-add_option(parser, c("-i", "--individual"), default=TRUE, action="store_true", help="Flag to output the classification for each individual")
parser<-add_option(parser, c("-m", "--metric"), default="qNet", help="the metric to use")
parser<-add_option(parser, c("-d", "--APpath"), default="../nanion8rand_AP_simulation/", help="Filepath to the metrics.rds")
args<-parse_args(parser)

set.seed(100)

#--- load libraries
library(rms)
library(ggplot2)
print(sessionInfo())

#--- get arguments
tdpfile<-args$tdpfile
useUQ<-args$uncertainty
omit<-args$omit
outputI<-args$individual
metric<-args$metric
APpath<-args$APpath

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
    try({ lrm(CiPA~get(metric), data=datadf, penalty=0, x=TRUE, y=TRUE, tol=tol, maxit=maxit) })
}

#--- read in drug TdP risk
drugtable<-read.csv(tdpfile)
drugnames<-as.character(drugtable$drug)

#--- do logistic regression
for(outdir in outdir_vec){
    # read in dataset
  if(metric =="qNet"){
    infile<-paste0(APpath,outdir,"metrics.rds")
    df<-readRDS(infile)
    df<-df[df$drug!="control" & df$dose!=0,]
    df<-df[!is.na(df$max_dv),] # depolarization failures
    df<-df[(df$dose)%in%c("1","2","3","4"),]                  #very important cause some drugs may have additional doses in metrics.rds!!!!
#print(head(df))
    if(input == "uncertainty"){
    qNettable<-aggregate(df[,metric], list(df$drug,df$sample), mean)  #name is qNettable, but it could be any metric
colnames(qNettable)<-c("drug","sample",metric)
     }else{
     qNettable<-aggregate(df[,metric], list(df$drug), mean)
colnames(qNettable)<-c("drug","qNet")
     qNettable$sample <- NA
     }
print(head(qNettable))
print(class(qNettable$qNet))
if(metric == "qNet")
qNettable$qNet<-qNettable$qNet/1000
    df<-merge(qNettable, drugtable[,c("drug","CiPA")], by.x="drug", by.y="drug")
   
 }
 
    #organize df
     df$drug<-factor(df$drug,levels=drugnames)
    df<-df[order(df$drug,df$sample),]
    df$CiPA<-ordered(df$CiPA)

    # perform logistic regression and leave-one-out cross validation for each drug
    errdf<-data.frame()
    probdf<-data.frame()
    allprobdf<-data.frame()   #corresponding to cvdf but without LOOCV
    cverrdf<-data.frame()
    cvprobdf<-data.frame()
    cvdf<-data.frame()                              #used to store individual sample's probs

        for(drug in c(NA,drugnames)){

            # fit logistic regression model
            if(is.na(drug)){
                print("training on all drugs")
                datadf<-df[, c("drug","CiPA","sample",metric)]
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
                    print("fitting all drugs failed? skipping...")
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
            y0<-testdf$CiPA
            if(is.na(drug)){
                probs<-predict(lmod, type="fitted.ind")
            }else{
                probs<-matrix(predict(lmod, newdata=testdf, type="fitted.ind"), ncol=length(levels(y0)))
            }
            
            #sometimes APD90 is NA, giving rise to NA in probs
            idx<- apply(probs,1,function(x) any(is.na(x)))
            probs[idx,]<- matrix(rep(c(0,0,1),sum(idx)),nrow=sum(idx),byrow=T) #this is to assume all NAs are due to repolarization failure, and thus high risk
            
            yPred<-apply(probs,1, function(x) which.max(x))
            pred_err<-mean(abs(yPred-as.integer(y0)))
          
            # append to data frame
            if(is.na(drug)){
                newrow<-data.frame( t(cfvec), error=pred_err)
                errdf<-rbind(errdf,newrow)
                
                allprobdf<-cbind(traindf[,1:3],probs); colnames(allprobdf)<-c("drug","CiPA","sample","low_prob","inter_prob","high_prob")
               
            }else{
                newrow<-data.frame(drug=drug,  t(cfvec), error=pred_err)
                cverrdf<-rbind(cverrdf,newrow)
               if(outputI) {                            #make detailed prob table here
                testdf<-cbind(testdf, high_prob=probs[,3], inter_prob=probs[,2], low_prob=probs[,1], pred_err=abs(yPred-as.integer(y0)))
                cvdf<-rbind(cvdf,testdf)
               }#if outputI
            }

            # detailed results
            newcols<-sapply(levels(y0), function(s) paste0("predict_",s), USE.NAMES=F)
            predictions<-factor(newcols[yPred], levels=newcols)
            tmpdf<-as.data.frame(tapply(1:nrow(testdf), list(drug=as.character(testdf$drug), predict=predictions), FUN=length))
            tmpdf[is.na(tmpdf)]<-0
            tmpdf$drug<-factor(rownames(tmpdf), levels=drugnames)
            
            tmpdf<-merge(drugtable[,c("drug","CiPA")],tmpdf)
            tmpdf<-tmpdf[order(tmpdf$drug),]
            rownames(tmpdf)<-NULL
            if(is.na(drug)){
                probdf<-rbind(probdf,tmpdf)
               
            }else{
                cvprobdf<-rbind(cvprobdf,tmpdf)
            }

        } # for drug
 

    outfile<-paste0(outdir,metric,"training_probs.csv")   #here training means "trained on all drugs"
    write.csv(probdf, outfile, row.names=F, quote=F)
    print(head(probdf))

    outfile<-paste0(outdir,metric,"training_errors.csv")
    write.csv(errdf, outfile, row.names=F, quote=F)
    print(head(errdf))

    outfile<-paste0(outdir,metric,"LOOCV_probs.csv")
    write.csv(cvprobdf, outfile, row.names=F, quote=F)
    print(head(cvprobdf))

    outfile<-paste0(outdir,metric,"LOOCV_errors.csv")
    write.csv(cverrdf, outfile, row.names=F, quote=F)
    print(head(cverrdf))
    
    if(outputI){
      outfile<-paste0(outdir,metric,"LOOCV_allprobes.csv")
    write.csv(cvdf, outfile, row.names=F, quote=F)
    print(head(cvdf))
      outfile<-paste0(outdir,metric,"training_allprobes.csv")
      write.csv(allprobdf,outfile,row.names=F,quote=F)
      }
}# for outdir
