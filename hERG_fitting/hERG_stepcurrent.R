proc_start<-proc.time()

options(warn=1)

library(optparse)

isWindows<-Sys.info()[["sysname"]]=="Windows"

#--- specify command line arguments
parser<-OptionParser()
parser<-add_option(parser, c("-d", "--drug"), default="all", type="character", help="Drug name [required]")
parser<-add_option(parser, c("-f", "--file"), default=NULL, type="character", help="Path to CSV file [required]")
parser<-add_option(parser, c("-p", "--pickedtime"), default=100, help="averaging the last p ms of the 10th sweep as the step current")

args<-parse_args(parser)

#--- required argument
if(is.null(args$drug)) stop("Missing drug argument!")
if(is.null(args$file)) stop("Missing file path argument!")

drug<-args$drug
drugfile<-args$file

# Verify file exists
if(!file.exists(drugfile)) stop(sprintf("File %s does not exist!", drugfile))

#--- optional argument
pickedtime<- args$pickedtime
#pickedtime<- 100


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

alldf<-extracteddf
write.table(alldf, "drug_block.csv", sep=",",quote=F,col.names=T,row.names=F)
