#Bnet calculation Github
hERG_par="/data/results_1s_HT_Bmax"
otherCh_path="/data/non_hERG/results"
cmaxfile<-"data/newCiPA1.csv"
DRUGNAMES=c("dofetilide","bepridil","sotalol","quinidine","cisapride","terfenadine","ondansetron","chlorpromazine","verapamil","ranolazine","mexiletine","diltiazem","disopyramide","ibutilide","vandetanib","azimilide","droperidol","domperidone","clarithromycin","risperidone","pimozide","clozapine","metoprolol","loratadine","tamoxifen","nifedipine","nitrendipine","astemizole")
Bnet2=vector()
for (drug in DRUGNAMES) {
params=block1=Ochparams=ICaLparams=INaLparams=INaparams=vector()
drugtable<-read.csv(cmaxfile)
all_herg_par<-read.csv(sprintf("%s/%s/IC50_samples.csv",hERG_par,drug))
otherCh=read.csv(sprintf("%s/%s/IC50_samples.csv",otherCh_path,drug))
print(sprintf("%s/%s/IC50_samples.csv",otherCh_path,drug))
cmax<-drugtable[tolower(as.character(drugtable$drug))==drug,"therapeutic"] # should be in nanomolar
doses=c(1,2,3,4)
for (i in 1:2000) {
	params[1]=all_herg_par[i,"hERG_IC50"]
	params[2]=all_herg_par[i,"hERG_h"]
	
	ICaLparams[1]=otherCh[i,"ICaL_IC50"]
	ICaLparams[2]=otherCh[i,"ICaL_h"]
	
	INaparams[1]=otherCh[i,"INa_IC50"]
	INaparams[2]=otherCh[i,"INa_h"]
	
	INaLparams[1]=otherCh[i,"INaL_IC50"]
	INaLparams[2]=otherCh[i,"INaL_h"]

        for(dose in doses){
            conc<-dose*cmax
			block=100*(1-1/(1+(conc/params[1])^params[2]))
			blockICaL=100*(1-1/(1+(conc/ICaLparams[1])^ICaLparams[2]))
			blockINa=100*(1-1/(1+(conc/INaparams[1])^INaparams[2]))
			blockINaL=100*(1-1/(1+(conc/INaLparams[1])^INaLparams[2]))
			Bnet=block-(blockICaL+blockINaL)
			block0=c(block,blockICaL,blockINa,blockINaL,Bnet,conc,dose,i,params,ICaLparams,INaparams,INaLparams)

block1=rbind(block1,block0)
		}
   }
   colnames(block1)=c("blockhERG","blockICaL","blockINa","blockINaL","Bnet","conc","dose","i","IhERG","hhERGparams","ICaLparams","hICaLparams","INaparams","hINaparams","INaLparams","hINaLparams")
   write.csv(block1,paste0("results/",drug,"_block2.csv"))
   Bnet1=c(drug,mean(block1[,5]),quantile(block1[,5], probs = c(0.025, 0.975), na.rm = FALSE, names = TRUE))
   Bnet2=rbind(Bnet2,Bnet1)
}
write.csv(Bnet2,paste0("results/","Bnet2.csv"))