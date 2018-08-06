# File:         keep_only_Li2017_IC50.R
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  R script to remove IC50 values from MCMC simulation outputs
#               that were not provided in the article
#               "Improving the in silico assessment of proarrhythmia risk by
#               combining hERG (Human Ether-à-go-go-Related Gene)
#               channel-drug binding kinetics and multichannel pharmacology"
#               by Zhihua Li, Sara Dutta, Jiansong Sheng, Phu N. Tran, Wendy
#               Wu, Kelly Chang, Thembi Mdluli, David G. Strauss, and Thomas
#               Colatsky
#               Link: http://circep.ahajournals.org/content/10/2/e004628.abstract
#

# read in Li et al. 2017 IC50 values
drugtable<-read.csv("data/Li2017_IC50.csv")

for(drug in drugtable$drug){
    drugdir<-sprintf("results/%s/",drug)
    print(drugdir)

    for(filename in c("IC50_optimal.csv", "IC50_samples.csv")){
        # read in IC50 table
        drugfile<-paste0(drugdir,filename)
        df<-read.csv(drugfile)

        # save original file for backup
        i<-0
        bkfile<-paste0(drugfile,".bkup.",i)
        while(file.exists(bkfile)){
            i<-i+1
            bkfile<-paste0(drugfile,".bkup.",i)
        }
        file.copy(drugfile, bkfile)

        # find columns with no IC50 value (set at 0) in Li et al. 2017
        pidx<-match(colnames(df), colnames(drugtable), nomatch=0)
        thistable<-drugtable[drugtable$drug==drug, pidx]
        isna<-which(is.na(thistable[1,]))

        # erase values and overwrite file with new table
        print(head(df[,isna]))
        df[,isna]<-NA
        write.csv(df, drugfile, row.names=F, quote=F)
    }#for filename
}#for drug
