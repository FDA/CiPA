# File:         stepprotocol.R
# Author:       Zhihua Li
# Date:         2016
# Version:      1.0
# 
# Description:  Helper R function to define voltage step protocols.
#

stepprotocol<-function(holdv, holdt,conv,cont,testv,testt,gapt){
    #generic step protocol consisting a holding step, a conditioning step, and a testing step
    #gapt is the gap from the end of the this pulse to the start of next pulse

    times<-c(0,holdt,cont+holdt,cont+holdt+testt)
    voltages<-c(holdv,conv,testv,holdv)          #note the unit is mv !!!
    eventdata<-data.frame(var="V",time=times,value=voltages,method="replace")

    peaktimes<-seq(holdt+cont,holdt+cont+testt,1)                                             #for obtaining the peak; should be within 1 s after depolarization
    finaltimes<-  holdt+cont+testt+gapt

    initialtime<- 0                               #included in times already
    outputtimes<-sort(unique(c(initialtime,peaktimes,finaltimes)))
    fulltimes<-sort(c(times, cleanEventTimes(outputtimes, times))) #times is actually eventtimes

    output<-list(fulltimes, peaktimes,eventdata)
}
