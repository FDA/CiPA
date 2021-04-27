# File:         get_pars.R
# Author:       Kelly Chang
# Date:         Sep 2017
# Version:      1.0
# 
# Description:  Helper R function that returns the model parameters for a
#               specified simulation type (e.g. control, drug with fixed
#               inputs, drug with uncertainty inputs, etc.).
#

get_pars<-function(parfile, drug="control", free=0, isamp=0, bootpath=NA, mcmcpath=NA, drop=NA, optimized=TRUE){
    if(free<0) stop("get_pars(): negative concentration!")
    if(isamp<0) stop("get_pars(): negative sample number!")
    if(drug!="control" && (is.na(bootpath) || is.na(mcmcpath)))
        stop("get_pars(): must specify bootpath and mcmcpath!")

    # read default parameters from file
    tmp<-read.table(parfile, col.names=c("param","value"))
    pars<-setNames(tmp$value, tmp$param)

    pars["GKrfc"] = 1.1                 #manual adjustment of IKr conductance
    pars["T"] = 37

    if(drug=="control"){
        pars["Kmax"]=0
        pars["Kt"]=0

        free=0
    }else{
        # load drug pars
        if(isamp==0){ # use optimal parameters
            hERGfile<- sprintf("%s/%s/pars.txt",bootpath,drug)
           if(file.exists(hERGfile)){
           tmp<-read.table(sprintf("%s/%s/pars.txt",bootpath,drug), header=FALSE, row.names=1)
             }else{                        #file not copied yet
             hERGpath<- paste0(bootpath,"/",drug)
             system(paste0("mkdir -p ",hERGpath))
             system(paste0("cp ~/documents/cardio/HERGMMCMAES/Output/",drug,"/noE4031subtraction/pars.txt ",hERGpath,"/"))
               tmp<-read.table(sprintf("%s/%s/pars.txt",bootpath,drug), header=FALSE, row.names=1)
            }
            rownames(tmp)[1]<-"Kmax"  #if file copied from HERGMMCMAES fitting, then the parameter is called "Kf": need to change!
            #and you need to change it here rather than within the brackets above!
            hERGpars<-t(tmp)[1,]
            
             tmp<-read.csv(sprintf("%s/%s/boot_pars.csv",bootpath,drug))
            tmp<-apply(tmp, 2, function(x) median(sample(x,2000,replace=T)))
             hERGpars<-unlist(tmp)
            
            IC50table<-read.csv(sprintf("%s/%s/IC50_optimal.csv",mcmcpath,drug))
        }else{ # use uncertainty sample
            tmp<-read.csv(sprintf("%s/%s/boot_pars.csv",bootpath,drug))
            tmp<-apply(tmp, 2, function(x) sample(x,2000,replace=T))
            hERGpars<-unlist(tmp[isamp,])
            mcmctable<-read.csv(sprintf("%s/%s/IC50_samples.csv",mcmcpath,drug))
            IC50table<-mcmctable[isamp,]
        }

        ispar<-names(hERGpars)%in%names(pars)
        pars[names(hERGpars)[ispar]]<-hERGpars[ispar]

        pars["Kt"]<-3.5e-5
        if(!is.na(drop)){
            if(drop=="hERG"){
                pars["Kmax"]=0
                pars["Kt"]=0
            }else{
                IC50table[,paste0(drop,"_IC50")]<-NA
                IC50table[,paste0(drop,"_h")]<-NA
            }
        }

        # IC50s
        if(free > 0 ){
            print("Scaling ionic currents according to IC50 values.")
            if("INa_IC50"%in%colnames(IC50table) && !is.na(IC50table$INa_IC50))
                pars["GNafc"]<- 1+(free/IC50table$INa_IC50)^IC50table$INa_h
            if("INaL_IC50"%in%colnames(IC50table) && !is.na(IC50table$INaL_IC50))
                pars["GNaLfc"]<-1+(free/IC50table$INaL_IC50)^IC50table$INaL_h
                
            if("ICaL_IC50"%in%colnames(IC50table) && !is.na(IC50table$ICaL_IC50))
                pars["PCafc"]<- 1+(free/IC50table$ICaL_IC50)^IC50table$ICaL_h
          
        }#if there is drug then do IC50s
    }

    if(optimized){
        print("Applying optimized scaling factors.")
        pars["GKrfc"]<-pars["GKrfc"]/1.114
        pars["GNaLfc"]<-pars["GNaLfc"]/2.661
        pars["GKsfc"]<-pars["GKsfc"]/1.870
        pars["GK1fc"]<-pars["GK1fc"]/1.698
        pars["PCafc"]<-pars["PCafc"]/1.007
    }

    # return new parameters
    pars
}
