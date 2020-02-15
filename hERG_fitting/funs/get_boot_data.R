# File:         get_boot_data.R
# Author:       Kelly Chang
# Date:         Sep 2017
# Version:      1.0
# 
# Description:  Helper R function which computes the mean fractional current
#               from Milnes protocol data. The data are resampled according
#               to the specified cell indices.
#

get_boot_data<-function(filepath, idx){
    # read original data
    datadf<-read.csv(filepath)
    expdf<-unique(datadf[,c("conc","exp")])
    expdf<-expdf[with(expdf,order(conc,exp)),]

    # bootstrap data
    cellsdf<-expdf[idx,]
    boot_list<-lapply(1:nrow(cellsdf),
                      function(i){
                          datadf[datadf$conc==cellsdf[i,"conc"] &
                                 datadf$exp==cellsdf[i,"exp"],
                                 c("time","sweep","frac")]}
                      )
    cellsdf$conc<-as.character(cellsdf$conc) # convert to characters for indexing

    # get mean traces with selected points
    fracdata<-list()
    for(conc in unique(cellsdf$conc)){
        expidx<-which(cellsdf$conc==conc)
        concdf<-do.call(rbind,boot_list[expidx])

        # compute mean of bootstrap data
        meandf<-aggregate(concdf[,"frac",drop=FALSE], by=list(time=concdf$time, sweep=concdf$sweep), FUN=mean)

        fracdata[[conc]]<-list()
        sweeps<-unique(meandf$sweep)
        if(any(sweeps!=1:length(sweeps))) stop("Bad sweeps!")
        for(i in sweeps){
            sweepdf<-meandf[meandf$sweep==i,c("time","frac")]
            sweepdf<-sweepdf[order(sweepdf$time),]

            # select starting point from first sweep
            if(i==1) startpt<-min(c(11,which(sweepdf$frac<=0.5)))
            fracdata[[conc]][[i]]<-sweepdf[startpt:nrow(sweepdf),]
        }# for each sweep
    }# for each dose

    fracdata
}
