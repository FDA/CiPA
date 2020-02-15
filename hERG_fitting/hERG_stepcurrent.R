proc_start<-proc.time()

options(warn=1)

library(optparse)

isWindows<-Sys.info()[["sysname"]]=="Windows"

#--- specify command line arguments
parser<-OptionParser()
parser<-add_option(parser, c("-d", "--drug"), default="all", type="character", help="Drug name [required]")
parser<-add_option(parser, c("-p", "--pickedtime"), default=100, help="averaging the last p ms of the 10th sweep as the step current")

args<-parse_args(parser)

#print(sessionInfo())

#--- required argument
if(is.null(args$drug)) stop("Missing drug argument!")
drug<-args$drug
drugvec<- c("dofetilide","bepridil","quinidine","sotalol","ibutilide", "azimilide", "disopyramide", "vandetanib",
             "cisapride","terfenadine","ondansetron","chlorpromazine","clarithromycin","risperidone","domperidone","astemizole", "pimozide","droperidol","clozapine",
             "ranolazine","mexiletine","diltiazem","verapamil","metoprolol","loratadine","tamoxifen","nifedipine","nitrendipine")

drugvec<-drug
#--- optional argument
pickedtime<- args$pickedtime
#pickedtime<- 100

alldf<-data.frame(drug=NA, conc=NA, units=NA, channel=NA, block=NA,exp=NA)
for(drug in drugvec){

drugfile<- paste0("data/",drug,".csv")
drugdf<-read.delim(drugfile,sep=",")
#only the last sweep is needed
drugdf<- drugdf[drugdf$sweep == 10,]

selectfun<-function(partdf){              #partdf is the chunk of df with a specific exp number and conc
   colnames(partdf)<-colnames(drugdf)
   lasttime<- which.max(partdf$time)
   includedpt <- partdf$time >= (lasttime - pickedtime)
   meanres<- mean(partdf$frac[includedpt])
   return(c(partdf$exp[1], partdf$conc[1], meanres))
   }
   
extractedlist<- by(drugdf, list(drugdf$exp, drugdf$conc), selectfun, simplify=T)
extracteddf<-do.call(rbind, extractedlist)
extracteddf<-as.data.frame(extracteddf)
colnames(extracteddf)<-c("exp","conc","block")
extracteddf$block <- 100*(1-extracteddf$block)
extracteddf$channel<- "hERG";  extracteddf$units<- "nM"; extracteddf$drug<- drug; 
alldf<-rbind(alldf,extracteddf)
}#drug vec

alldf<-alldf[-1,]
write.table(alldf, "drug_block.csv", sep=",",quote=F,col.names=T,row.names=F)

