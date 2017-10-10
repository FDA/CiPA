# File:         process_patchclamp_data.R
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  R script for processing patch clamp data from Li et al. 2017
#               to be used for Markov-chain Monte Carlo simulation and
#               uncertainty quantification.
#               This data was used in the study "Improving the in silico
#               assessment of proarrhythmia risk by combining hERG (Human
#               Ether-à-go-go-Related Gene) channel-drug binding kinetics
#               and multichannel pharmacology"
#               by Zhihua Li, Sara Dutta, Jiansong Sheng, Phu N. Tran, Wendy
#               Wu, Kelly Chang, Thembi Mdluli, David G. Strauss, and Thomas
#               Colatsky
#               Link: http://circep.ahajournals.org/content/10/2/e004628.abstract
#


# read in patch clamp data used for Li et al. 2017
datadf<-read.csv("mergedpatchclampdata-20160514.csv")

# convert to lowercase
colnames(datadf)<-tolower(colnames(datadf))
datadf$drug<-tolower(datadf$drug)

# remove Crumb et al. 2016 hERG data and add White Oak data
datadf<-datadf[datadf$channel!="hERG",]
tmpdf<-read.csv("WO_hERG_block.csv")
colnames(tmpdf)<-tolower(colnames(tmpdf))
tmpdf$drug<-tolower(tmpdf$drug)
tmpdf$pacing<-NA
datadf<-rbind(datadf,tmpdf)

# convert micromolar to nanomolar
isMicro<-datadf$Units=="uM"
datadf[isMicro,"conc"]<-datadf[isMicro,"conc"]*1e3
datadf[isMicro,"units"]<-"nM"

# convert channel names to current names
datadf$channel<-sapply(datadf$channel, function(x) gsub(" ","_",x))
channels<-as.character(datadf$channel)
channels[channels=="Calcium"]<-"ICaL"
channels[channels=="Kv4.3"]<-"Ito"
channels[channels=="Peak_sodium"]<-"INa"
channels[channels=="Late_sodium"]<-"INaL"
datadf$channel<-channels

# save to new file
write.csv(datadf, "drug_block.csv", row.names=F, quote=F)
