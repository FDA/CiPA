# File:         compute_qNet_CI.R
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  R script to compute credible intervals for qNet and plot
#               the results for fixed-input and uncertainty-input simulations.
#               For help, run this script with command line option "-h".
#

#--- specify command line arguments
library(optparse)
parser<-OptionParser()
parser<-add_option(parser, c("-t", "--tdpfile"), default="data/CiPA_training_drugs.csv", help="Filepath to table of known TdP risk levels [default data/CiPA_training_drugs.csv]")
parser<-add_option(parser, c("-a", "--alpha"), default=0.05, help="Significance level [default 0.05]")

args<-parse_args(parser)

#--- load libraries
library(ggplot2)

#--- get arguments
tdpfile<-args$tdpfile
alpha<-args$alpha

# read results
outdir<-"results/"
infile<-paste0(outdir,"fixed/metrics.rds")
fixdf<-readRDS(infile)
infile<-paste0(outdir,"uncertainty/metrics.rds")
sampdf<-readRDS(infile)

# compute quantiles
fixdf<-fixdf[fixdf$drug!="control",]
cidf<-data.frame()
for(drug in unique(fixdf$drug)){
    for(dose in sort(unique(fixdf$dose))){
        print(sprintf("%s: %gx Cmax",drug,dose))

        # fixed value
        thisfixdf<-fixdf[fixdf$drug==drug & fixdf$dose==dose & is.na(fixdf$sample),]
        fix_val<-ifelse(is.na(thisfixdf$max_dv), NA, thisfixdf$qNet) # depolarization failure

        # sample intervals
        thisdf<-sampdf[sampdf$drug==drug & sampdf$dose==dose,]
        thisdf<-thisdf[order(thisdf$sample),] # order samples
        if(any(thisdf$sample!=1:nrow(thisdf)))
            stop("Missing samples!")
        rownames(thisdf)<-NULL

        # samples with depolarization failures
        isnoDepol<-is.na(thisdf$max_dv)
        thisdf[isnoDepol,"qNet"]<-NA

        # compute quantiles and save to dataframe
        ci_lims<-quantile(thisdf$qNet,probs=c(alpha/2,1-alpha/2),na.rm=TRUE)
        out<-data.frame(drug=drug, dose=dose, value=fix_val, samples=nrow(thisdf), conf=1-alpha, lower=ci_lims[[1]], upper=ci_lims[[2]])
        cidf<-rbind(cidf,out)
    } # for dose
} # for drug
print(head(cidf))

outfile<-(sprintf("results/qNet_CI%g.csv",(1-alpha)*100))
write.csv(cidf, outfile, row.names=F, quote=F)

# prepare results for plotting
conf<-1-alpha
cidf[,c("value","lower","upper")]<-cidf[,c("value","lower","upper")]/1e3 # convert units

# get risk categories
drugtable<-read.csv(tdpfile)
cidf<-merge(cidf, drugtable[,c("drug","CiPA")], all.x=T)
cidf$CiPA<-factor(cidf$CiPA, levels=c(2,1,0))

# data frame of drug labels
newdf<-aggregate(dose ~ drug, data=cidf, max)
colnames(newdf)<-c("drug","max")
newdf<-merge(cidf,newdf,by.x="drug",by.y="drug",all=TRUE)
newdf<-newdf[newdf$dose==newdf$max,]
tofix<-which(is.na(newdf$value))
for(fixi in tofix)
    newdf[fixi,"value"]<-mean(as.numeric(newdf[fixi, c("lower","upper")]))

# plot results
system(paste0("mkdir -p figs"))
figfile<-sprintf("figs/qNet_CI%g.pdf",conf*100)
pdf(figfile, width=8, height=4.5)

p<-ggplot(cidf, aes(x=dose, y=value, group=drug))
p<-p+scale_color_brewer(palette="Set1")
p<-p+scale_fill_brewer(palette="Set1")
p<-p+geom_line(aes(color=CiPA))
p<-p+geom_ribbon(aes(ymin=lower, ymax=upper, fill=CiPA, group=drug), alpha=0.3)
p<-p+geom_text(data=newdf, aes(label=drug, color=CiPA), hjust=-0.2, vjust=0.5, size=3, show.legend=F)
p<-p+coord_cartesian(xlim=c(0,max(cidf$dose)*1.2))
p<-p+xlab(~"Concentration (\u00D7"~C[max]*")")
p<-p+ylab(~"qNet ("*mu*"C/"*mu*"F)")
p<-p+theme_bw()
p<-p+theme(legend.position="none")
print(p)
dev.off()
