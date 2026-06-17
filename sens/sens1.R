rm(list=ls())


#mei 28 drug t1 t2
library(rms)
set.seed(1)
library(ggplot2)

scale<-"sara"
CL<-2000
metric<-"qNet"
dataset<-"chantest"

last_dose<-4
dose<-last_dose
outdir<-sprintf("results/%s/scale_%s/CL%d/",dataset,scale,CL)
last_dose<-4
fitmethod<-"lrm"
#fitmethod<-"ordLORgee"
#fitmethod<-"repolr"

colorvec<-c("red","dark green","blue","black","cyan","purple","brown")
colordf<-data.frame(
   drug = c("dofetilide","bepridil","quinidine","sotalol","ibutilide", "azimilide", "vandetanib",
             "cisapride","terfenadine","ondansetron","chlorpromazine","clarithromycin","risperidone","astemizole", "pimozide","droperidol","clozapine",
             "ranolazine","mexiletine","diltiazem","verapamil","metoprolol","loratadine","tamoxifen","nifedipine","nitrendipine"),
   classidx=c(rep(2,7),rep(1,10),rep(0,9)),
   coloridx=c(rep(1,7),rep(3,10),rep(2,9)),
   isTraining=c(1,1,1,1,0,0,0,
                1,1,1,1,0,0,0,0,0,0,
                1,1,1,1,0,0,0,0,0)
 )


system(paste0("mkdir -p ",outdir))
figdir<-sub("results","figs",outdir)
system(paste0("mkdir -p ",figdir))

source("fit_funs.R")
if(fitmethod=="polr"){
    library(MASS)
    try_fit<-try_polr
}else if(fitmethod=="lrm"){
    library(rms)
    try_fit<-try_lrm
}else if(fitmethod=="repolr"){
    require(repolr)
    try_fit<-try_repolr
}else if(fitmethod=="ordLORgee"){
    require(multgee)
    try_fit<-try_ordLORgee
}

# read in dataset

thrdf<-data.frame()

if(metric !="model5"){
df<- as.data.frame(readRDS(paste0("../",dataset,"_AP_simulation/results/uncertainty/metrics.rds")))  
#mei head
#df<-df[df$drug!="disopyramide"&df$drug!="domperidone",]
#mei over
df<-df[df$drug!="control" & df$dose>0,]
todrop<-is.na(df$max_dv)#df$noDepols==250
print(sprintf("Number of depolarization failures: %d/%d",sum(todrop),nrow(df)))
df<-df[!todrop,]
df$class<-df$max_dv>0
df$qNet <- df$qNet/1000
print(sprintf("Number of EADs: %d/%d",sum(df$class),nrow(df)))

if(dataset=="original"){
    nboots<-0
}else{
    #nboots<-max(df$boot) # using all bootstraps
    nboots<-2000
}

cols<-c("drug","sample")


# create thresholds for each range of doses
print(sprintf("%s, %s scale, CL = %d ms, nboot = %d, dose <= %d Cmax",dataset,scale,CL,nboots,dose))

brows<-df$dose<=dose
wtdf<-aggregate(df[brows,metric,drop=FALSE], by=as.list(df[brows,cols]), mean) # equally weighted across doses
wtdf$drug<-as.factor(wtdf$drug)
}else{                            #if metric == model5
   if(dataset=="pureHTS"){                   #shouldn't use "hybrid" as dataset here!
    Hilldataset<- "pureHTS"
   }else if(dataset=="pureManual"){
    Hilldataset<-"pureManual"          #note I didn't use validationManual_Hill_fitting here
   }
   wtdf<-data.frame(drug=rep(colordf$drug, each=2000),sample=rep(1:2000,28),model5=NA)
   for(drug in colordf$drug){
      IC50table<- read.delim(paste0("/scratch/lizhi/cardio/UQ/",Hilldataset,"_Hill_fitting/results/",drug,"/IC50_samples.csv"),sep=",")
      wtdf$model5[wtdf$drug==drug]<- with(IC50table, -log(ICaL_IC50) - -log(hERG_IC50))
    }#for drug




}#if metric !=model5
# get 95% CI
q2.5<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
                                                  FUN=function(x) quantile(x,probs=0.025,na.rm=TRUE))
q97.5<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
                                                  FUN=function(x) quantile(x,probs=0.975,na.rm=TRUE))
                                                  
q50<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
                                                  FUN=function(x) quantile(x,probs=0.5,na.rm=TRUE))                                                  
cidf<-merge(q2.5,q97.5,by=c("drug"),suffixes=c("_0.025","_0.975"))
cidf<-merge(cidf, q50, by="drug")
colnames(cidf)[4]<-paste0(metric,"_0.5")
outdf<-cidf

outfile<-sprintf("%s/%s_allCIthresholds.rds",outdir,metric)
saveRDS(outdf, outfile)


outdf<- merge(colordf, outdf,by="drug")
cidf<-outdf

desiredorder<-drugnames<-c("ibutilide","vandetanib","bepridil","azimilide","dofetilide","quinidine","sotalol",  "disopyramide", 
             "domperidone","pimozide","droperidol","cisapride","terfenadine","astemizole","ondansetron","clarithromycin","chlorpromazine","clozapine","risperidone", 
             "tamoxifen","loratadine","verapamil","metoprolol","nitrendipine","diltiazem","ranolazine","nifedipine","mexiletine")

outdf$drug<- factor(outdf$drug, levels=rev(desiredorder))
outdf$class<- 2-outdf$classidx
lower<-paste(metric,0.025,sep="_")
upper<-paste(metric,0.975,sep="_")
middle<-paste(metric,0.5,sep="_")
lower_bds<-aggregate(outdf[,lower], by=list(risk=outdf$class), FUN=min)
upper_bds<-aggregate(outdf[,upper], by=list(risk=outdf$class), FUN=max)
lower_bds<-lower_bds[match(c(0:2),lower_bds$risk),]
upper_bds<-upper_bds[match(c(0:2),upper_bds$risk),]
thresholds<-apply(cbind(head(upper_bds$x,-1),tail(lower_bds$x,-1)), 1, max)
thresholds<-sort(thresholds)
names(thresholds)<-c("High","Intermediate")

newrow<-data.frame(last_dose=dose, type="95CI", t(thresholds))
thrdf<-rbind(thrdf,newrow)

figfile<-sub("results","figs",sub(".rds",".pdf",outfile))
pdf(figfile, width=8, height=4)
p<-ggplot(outdf, aes(x=drug, color=as.character(class)))    #ggplot2 doesn't like color idx being numbers
p<-p+geom_errorbar(aes_string(ymin=lower, ymax=upper))
for(thresh in thresholds)
    p<-p+geom_hline(yintercept=thresh, linetype="dotted")
p<-p+ylab(sprintf("%s_score_1-%gX_Cmax",metric,dose))
p<-p+coord_flip()
p<-p+scale_color_brewer(NULL, palette="Set1")
p<-p+theme_bw()
print(p)
dev.off()

# get ordinal logistic regression thresholds
wtdf<-merge(wtdf,colordf)

wtdf$class <- ordered(2-wtdf$classidx)  #change the order because in Kelly's code "high" is first class (0).
                                        #but in mycompute_TdP_error.R "high" is 2 as planned

# fit all data using single predictor
datadf<-wtdf[,c("drug","class","sample",metric)]
#datadf$drug<- factor(datadf$drug, levels=rev(desiredorder))
maxsamp <- 2000
datadf<-datadf[datadf$sample<=maxsamp,]
maxreplacedidx<-0        #to break the perfect correlation within drugs for GEE methods
                                               #for ordinary lrm can use 0
if(maxreplacedidx !=0){                                               
test<-do.call(rbind,by(datadf, datadf$drug, function(x) {
                                      class<-x$class; fullvec<-0:2;
                                         for(i in 1:maxreplacedidx){
                                             c<-class[i]; idx<-fullvec%in%c;class[i]<-sample(fullvec[!idx],1);
                                             }
                                        x$class<-class;return(x)
                                             }))
datadf<-test[order(test$drug,test$sample),]  #sort by sample (time) is required by repolr but not ordLORgee
}


if(fitmethod=="lrm"){
dd<-datadist(datadf)
options(datadist="dd")
}else if(fitmethod=="repolr"){             #repolr doesn't like factors! and doesn't like class being 0?

datadf$class<-as.integer(datadf$class)  #if -1 then the values are the same as factor (0,1,2)
}
lmod<-try_fit(datadf, metric, 0) # no penalty
print(lmod)

if(inherits(lmod, "try-error")){
    print(sprintf("fitting dose %d failed! skipping...",dose))
    next
}

# save coefficients
if(fitmethod=="polr"){
    print(sprintf("Convergence code: %d, Number of iterations: %d",lmod$convergence,lmod$niter))
    cf<-lmod$coefficients[[metric]]
    ints<-lmod$zeta
    cfvec<-c()
    for(kint in 1:length(ints))
        cfvec[[paste0("intercept",kint)]]<-ints[[kint]]
    cfvec[["slope"]]<-cf
}else{
    print(sprintf("Convergence failure: %s",lmod$fail))
    cf<-coefficients(lmod)
    cfvec<-c()
    for(kint in 1:(length(cf)-1))
        cfvec[[paste0("intercept",kint)]]<-cf[[kint]]
    cfvec[["slope"]]<-cf[[length(cf)]]
}
t1<- -cfvec[["intercept1"]]/cfvec[["slope"]]
t2<- -cfvec[["intercept2"]]/cfvec[["slope"]]
#use formal math
ebeta1<-exp(-cfvec[["intercept1"]]); ebeta2<- exp(-cfvec[["intercept2"]])
t1 <- log(ebeta1*ebeta2/-(2*ebeta1-ebeta2))/cfvec[["slope"]]
t2 <- log(exp(-cfvec[["intercept2"]])-2*exp(-cfvec[["intercept1"]]))/cfvec[["slope"]]

thresholds<-c(High=t2, Intermediate=t1)
print(thresholds)

high1<-t2;intermediate1<-t1;

#meih


#meio
#mei 28 drug t1 t2 for 
#meih

xmdf<- as.data.frame(readRDS(paste0("../",dataset,"_AP_simulation/results/uncertainty/metrics.rds")))

xmthreshold1<-list()#28drug
xmthreshold2<-list()
xmdrug<-unique(xmdf$drug)

for (xmn in 1:length(xmdrug)){
thresholdslist1<-list()#16thresholds
thresholdslist2<-list()
xmnum<-seq(0.9,1.1,by=0.1)

for (m in 1:length(xmnum)){

xm1<-xmdf[xmdf$drug!=xmdrug[xmn],]
xm2<-xmdf[xmdf$drug==xmdrug[xmn],]
xm2$qNet<-xm2$qNet*xmnum[m]
#meio
xmwtdf<-rbind(xm1,xm2)


source("fit_funs.R")
if(fitmethod=="polr"){
    library(MASS)
    try_fit<-try_polr
}else if(fitmethod=="lrm"){
    library(rms)
    try_fit<-try_lrm
}else if(fitmethod=="repolr"){
    require(repolr)
    try_fit<-try_repolr
}else if(fitmethod=="ordLORgee"){
    require(multgee)
    try_fit<-try_ordLORgee
}

# read in dataset

thrdf<-data.frame()

if(metric !="model5"){
df<-xmwtdf  
df<-df[df$drug!="control" & df$dose>0,]
todrop<-is.na(df$max_dv)#df$noDepols==250
print(sprintf("Number of depolarization failures: %d/%d",sum(todrop),nrow(df)))
df<-df[!todrop,]
df$class<-df$max_dv>0
df$qNet <- df$qNet/1000
print(sprintf("Number of EADs: %d/%d",sum(df$class),nrow(df)))


if(dataset=="original"){
    nboots<-0
}else{
    #nboots<-max(df$boot) # using all bootstraps
    nboots<-2000
}

cols<-c("drug","sample")


# create thresholds for each range of doses
print(sprintf("%s, %s scale, CL = %d ms, nboot = %d, dose <= %d Cmax",dataset,scale,CL,nboots,dose))

brows<-df$dose<=dose
wtdf<-aggregate(df[brows,metric,drop=FALSE], by=as.list(df[brows,cols]), mean) # equally weighted across doses
wtdf$drug<-as.factor(wtdf$drug)
}else{                            #if metric == model5
   if(dataset=="pureHTS"){                   #shouldn't use "hybrid" as dataset here!
    Hilldataset<- "pureHTS"
   }else if(dataset=="pureManual"){
    Hilldataset<-"pureManual"          #note I didn't use validationManual_Hill_fitting here
   }
   wtdf<-data.frame(drug=rep(colordf$drug, each=2000),sample=rep(1:2000,28),model5=NA)
   for(drug in colordf$drug){
      IC50table<- read.delim(paste0("/scratch/lizhi/cardio/UQ/",Hilldataset,"_Hill_fitting/results/",drug,"/IC50_samples.csv"),sep=",")
      wtdf$model5[wtdf$drug==drug]<- with(IC50table, -log(ICaL_IC50) - -log(hERG_IC50))
    }#for drug




}#if metric !=model5
# get 95% CI
q2.5<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
                                                  FUN=function(x) quantile(x,probs=0.025,na.rm=TRUE))
q97.5<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
                                                  FUN=function(x) quantile(x,probs=0.975,na.rm=TRUE))
                                                  
q50<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
                                                  FUN=function(x) quantile(x,probs=0.5,na.rm=TRUE))                                                  
cidf<-merge(q2.5,q97.5,by=c("drug"),suffixes=c("_0.025","_0.975"))
cidf<-merge(cidf, q50, by="drug")
colnames(cidf)[4]<-paste0(metric,"_0.5")
outdf<-cidf

outfile<-sprintf("%s/%s_allCIthresholds.rds",outdir,metric)
saveRDS(outdf, outfile)


outdf<- merge(colordf, outdf,by="drug")
cidf<-outdf

desiredorder<-drugnames<-c("ibutilide","vandetanib","bepridil","azimilide","dofetilide","quinidine","sotalol",  "disopyramide", 
             "domperidone","pimozide","droperidol","cisapride","terfenadine","astemizole","ondansetron","clarithromycin","chlorpromazine","clozapine","risperidone", 
             "tamoxifen","loratadine","verapamil","metoprolol","nitrendipine","diltiazem","ranolazine","nifedipine","mexiletine")


outdf$drug<- factor(outdf$drug, levels=rev(desiredorder))
outdf$class<- 2-outdf$classidx
lower<-paste(metric,0.025,sep="_")
upper<-paste(metric,0.975,sep="_")
middle<-paste(metric,0.5,sep="_")
lower_bds<-aggregate(outdf[,lower], by=list(risk=outdf$class), FUN=min)
upper_bds<-aggregate(outdf[,upper], by=list(risk=outdf$class), FUN=max)
lower_bds<-lower_bds[match(c(0:2),lower_bds$risk),]
upper_bds<-upper_bds[match(c(0:2),upper_bds$risk),]
thresholds<-apply(cbind(head(upper_bds$x,-1),tail(lower_bds$x,-1)), 1, max)
thresholds<-sort(thresholds)
names(thresholds)<-c("High","Intermediate")

newrow<-data.frame(last_dose=dose, type="95CI", t(thresholds))
thrdf<-rbind(thrdf,newrow)

figfile<-sub("results","figs",sub(".rds",".pdf",outfile))
pdf(figfile, width=8, height=4)
p<-ggplot(outdf, aes(x=drug, color=as.character(class)))    #ggplot2 doesn't like color idx being numbers
p<-p+geom_errorbar(aes_string(ymin=lower, ymax=upper))
for(thresh in thresholds)
    p<-p+geom_hline(yintercept=thresh, linetype="dotted")
p<-p+ylab(sprintf("%s_score_1-%gX_Cmax",metric,dose))
p<-p+coord_flip()
p<-p+scale_color_brewer(NULL, palette="Set1")
p<-p+theme_bw()
print(p)
dev.off()

# get ordinal logistic regression thresholds
wtdf<-merge(wtdf,colordf)

wtdf$class <- ordered(2-wtdf$classidx)  #change the order because in Kelly's code "high" is first class (0).
                                        #but in mycompute_TdP_error.R "high" is 2 as planned

# fit all data using single predictor
datadf<-wtdf[,c("drug","class","sample",metric)]
#datadf$drug<- factor(datadf$drug, levels=rev(desiredorder))
maxsamp <- 2000
datadf<-datadf[datadf$sample<=maxsamp,]
maxreplacedidx<-0                           #to break the perfect correlation within drugs for GEE methods
                                               #for ordinary lrm can use 0
if(maxreplacedidx !=0){                                               
test<-do.call(rbind,by(datadf, datadf$drug, function(x) {
                                      class<-x$class; fullvec<-0:2;
                                         for(i in 1:maxreplacedidx){
                                             c<-class[i]; idx<-fullvec%in%c;class[i]<-sample(fullvec[!idx],1);
                                             }
                                        x$class<-class;return(x)
                                             }))
datadf<-test[order(test$drug,test$sample),]  #sort by sample (time) is required by repolr but not ordLORgee
}


if(fitmethod=="lrm"){
dd<-datadist(datadf)
options(datadist="dd")
}else if(fitmethod=="repolr"){             #repolr doesn't like factors! and doesn't like class being 0?

datadf$class<-as.integer(datadf$class)  #if -1 then the values are the same as factor (0,1,2)
}
lmod<-try_fit(datadf, metric, 0) # no penalty
print(lmod)

if(inherits(lmod, "try-error")){
    print(sprintf("fitting dose %d failed! skipping...",dose))
    next
}

# save coefficients
if(fitmethod=="polr"){
    print(sprintf("Convergence code: %d, Number of iterations: %d",lmod$convergence,lmod$niter))
    cf<-lmod$coefficients[[metric]]
    ints<-lmod$zeta
    cfvec<-c()
    for(kint in 1:length(ints))
        cfvec[[paste0("intercept",kint)]]<-ints[[kint]]
    cfvec[["slope"]]<-cf
}else{
    print(sprintf("Convergence failure: %s",lmod$fail))
    cf<-coefficients(lmod)
    cfvec<-c()
    for(kint in 1:(length(cf)-1))
        cfvec[[paste0("intercept",kint)]]<-cf[[kint]]
    cfvec[["slope"]]<-cf[[length(cf)]]
}
t1<- -cfvec[["intercept1"]]/cfvec[["slope"]]
t2<- -cfvec[["intercept2"]]/cfvec[["slope"]]
#use formal math
ebeta1<-exp(-cfvec[["intercept1"]]); ebeta2<- exp(-cfvec[["intercept2"]])
t1 <- log(ebeta1*ebeta2/-(2*ebeta1-ebeta2))/cfvec[["slope"]]
t2 <- log(exp(-cfvec[["intercept2"]])-2*exp(-cfvec[["intercept1"]]))/cfvec[["slope"]]

#thresholds<-c(High=t2, Intermediate=t1)

#meih



thresholdslist1[[m]]<-c(xmdrug[xmn],xmnum[m],t1,intermediate1,"threshold2")

thresholdslist2[[m]]<-c(xmdrug[xmn],xmnum[m],t2,high1,"threshold1")
print(thresholdslist1[[m]])

}
xmthreshold1[[xmn]]<-do.call(rbind,thresholdslist1)
xmthreshold2[[xmn]]<-do.call(rbind,thresholdslist2)
}

xmthreshold11<-do.call(rbind,xmthreshold1)
xmthreshold12<-do.call(rbind,xmthreshold2)
xmthreshold13<-rbind(xmthreshold11,xmthreshold12)


print(xmthreshold13)
print(xmdrug)
print(intermediate1)
print(high1)
write.csv(xmthreshold13,"28_sens.csv")
