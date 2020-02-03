rm(list=ls())

library(reshape2)

drugnames<-c("dofetilide","bepridil","sotalol","quinidine","cisapride","terfenadine","ondansetron","chlorpromazine","verapamil","ranolazine","mexiletine","diltiazem")
conf<-0.95 # confidence level

df<-data.frame()
for(drug in drugnames){
    print(drug)
    outdir<-sprintf("results/%s/",drug)
    tmpdf<-read.csv(sprintf("%s/boot_CI%g.csv",outdir,conf*100))
    tmpdf$drug<-drug

    # make long table
    tmpdf2<-melt(tmpdf, id.vars=c("drug","param"))
    tmpdf2$label<-mapply(function(str1,str2) paste0(str1,"_",str2), tmpdf2$param, tmpdf2$variable)

    # make wide table
    tmpdf3<-dcast(tmpdf2[,c("drug","label","value")], drug+conf~label, FUN=mean)
    df<-rbind(df,tmpdf3)
} # for drug

params<-c("log10Kmax","log10Ku","n","log10halfmax","Vhalf","log10slope")
variables<-c("value","lower","upper")
coldf<-expand.grid(variable=variables, param=params)
cols<-mapply(function(str1,str2) paste0(str1,"_",str2), coldf$param, coldf$variable)
df<-df[,c("drug",cols)]
write.csv(df, sprintf("results/boot_CI%g_summary.csv",conf*100),row.names=F)
