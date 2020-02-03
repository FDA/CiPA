source("compute_likelihood8.R")
source("calc_lsq_ig4.R")
source("makeMasterFrameFaster.R")
source("calculate_real_acceptance_rate.R")
library(FME)
library(optparse)

#--- specify command line arguments
parser<-OptionParser()
parser<-add_option(parser, c("-d", "--drug"), type="character", help="Drug name [required]")
#parser<-add_option(parser, c("-f", "--filepath"), default="data/drug_block.csv", help="Path to file containing patch clamp fractional block data for the specified drug [default data/drug_block.csv].")
parser<-add_option(parser, c("-s", "--seed"), default=100, type="integer", help="Random seed [default 100]")
parser<-add_option(parser, c("-n", "--nsamples"), default=2000, type="integer", help="Number of samples to collect [default 2000]")
parser<-add_option(parser, c("-b", "--burnin"), default=10000, type="integer", help="Burnin length [default 10000]")
#parser<-add_option(parser, c("-b", "--burnin"), default=10000, type="integer", help="Burnin length [default 10000]")

parser<-add_option(parser, c("-t", "--thin"), default=10, type="integer", help="Thinning rate [default 10]")
parser<-add_option(parser, c("-c", "--channel"), type="character", default="all", help="Channel or Current name")
args<-parse_args(parser)
ti=1
#datafile<-args$filepath
seednum <- args$seed
nsamp <- args$nsamples
burnin <- args$burnin
thin<-args$thin
targetDrug <- args$drug
targetChannel <- args$channel
set.seed(seednum)
targetChannel="hERG"   
#targetDrug <- "risperidone"  #iman

print("Loading data.")
print(sprintf("Random seed = %d.",seednum))
#iman .................................................here in make... he is making the initial csv file ...
myframe <- makeMasterFrameFaster(trim = 1.0,
		#siteExclusionList = c("Nanion","Bayer","Crumb","FDA"),
#		siteExclusionList = c("NanionJapanPhysiologicalSyncropatch",
#				"NanionMunichPhysiologicalPatchliner",
#				"NanionPatchlinerPhysiological",
#				"CRLChantestMilnes",
#				"Bayer",
#				"Crumbmanual",
#				"FDAmanual"),
		siteExclusionList = NULL,
		targetChannel,
		targetDrug)
#print(head(myframe))
print("Calculating initial guess.")
#iman ----------------------------------------------------- Initial guest calculatio by lsq -------
startingVector <- calc_lsq_ig4(myframe) 
#iman output is 5 parameters +h&I for each site  #based on Dr li advice , compare initials with gary
#--------------------------------------------------------------------------------------------------



print("Success in calculating initial guess.")

myControl <- list()
myControl$factr <- 1e9
myControl$maxit <- 1000


#iman
Itr_UL=500000
Total_itr=0
#Generate lower and upper bounds.
print("Generating lower and upper bounds.")
g <- names(startingVector)
LB <- rep(-Inf,length(g))
UB <- rep(Inf,length(g))
names(LB) <- g
names(UB) <- g
getHill <- grepl("_Hill$",g)
LB[getHill] <- 0.5
UB[getHill] <- 2.0
LB["globalsigma"] <-  0.001
UB["globalsigma"] <- 50.0

pIC50_scale_tag <- grepl("_pIC50_scale$",g)
LB[pIC50_scale_tag] <- 0.001
#LB[pIC50_scale_tag] <- -Inf
UB[pIC50_scale_tag] <- Inf

pIC50_location_tag <- grepl("_pIC50_location$",g)
LB[pIC50_location_tag] <- -3.99   #iman why is brad using -4 while his unit is molar!
UB[pIC50_location_tag] <- Inf

Hill_shape_tag <- grepl("_Hill_shape$",g)
LB[Hill_shape_tag] <- 0.001
UB[Hill_shape_tag] <- Inf

Hill_scale_tag <- grepl("_Hill_scale$",g)
LB[Hill_scale_tag] <- 0.001
UB[Hill_scale_tag] <- Inf


#Load the top-level hierarchy gamma parameters.
#iman here mu value is 8.126953 while it seems it should be around 12 , ? like starting vector 
priorParameters <- read.csv("../../Input/Gamma_Top_Level/Gamma_Top_Level.csv") 
mySolution <- optim(par = startingVector,
		fn = compute_likelihood8, #iman inside the compute_likelihood8, he used -4 again
		inputFrame = myframe,
		negation = -2, #iman it is -1 inside the fn
		printData = FALSE,
		pIC50LocationShift = -4, #iamn it has been defined in side of fn , is it neccessary to redefine it?
		priorParameters = priorParameters,
		method = "L-BFGS-B",
		control = myControl,
		lower = LB,
		upper = UB,
		hessian = TRUE)

comparisonTibble <- data.frame(fieldnames = names(startingVector),
		ig = startingVector,
		MLE = mySolution$par)

uniqueSites <- unique(myframe$sitename)
myDrug <- unique(myframe$drug)
mleVector <- mySolution$par

print("startingVector")
print(startingVector)

print("mleVector")
print(mleVector)

#pdfFilename <- sprintf("%s-%s-%s-test.pdf",targetDrug,targetChannel,seednum)
#pdf(pdfFilename)
#for(i in 1:length(uniqueSites))
#{
#	
#	
#	subframe <- subset(myframe,sitename == uniqueSites[i])
#	scaledConcs <- subframe$conc / 1e9
#	#iman
#	scaledConcs=scaledConcs*10^6
#	xlower <- 10^(floor(log10(min(scaledConcs))) - 1)
#	xupper <- 10^(ceiling(log10(max(scaledConcs))) + 1)
#	
#	plot(scaledConcs,
#			subframe$block,
#			log = "x",
#			main = sprintf("%s-%s-%s",uniqueSites[i],unique(subframe$drug),targetChannel),
#			xlim = c(xlower,xupper),
#			ylim = c(-50,150),
#			xlab = "Concentration, x [molar]",
#			ylab = "Percentage Block")
#	
#	xSmooth <- 10^seq(log10(xlower),log10(xupper),length.out = 500)
#	
#	#Least-squares parameters extraction from startingVector.
#	pIC50_name_tag <- paste0(myDrug,"_",targetChannel,"_",uniqueSites[i],"_pIC50")
#	Hill_name_tag <- paste0(myDrug,"_",targetChannel,"_",uniqueSites[i],"_Hill")
#	
#	pIC50.lsq <- startingVector[pIC50_name_tag]
#	IC50.lsq <- 10^(6 - pIC50.lsq)
#	Hill.lsq <- startingVector[Hill_name_tag]
#	
#	pIC50.MLE <- mleVector[pIC50_name_tag]
#	IC50.MLE <- 10^(6 - pIC50.MLE)
#	Hill.MLE <- mleVector[Hill_name_tag]
#	
#	
#	ySmooth.leastsquares <- 100 / (1 + (IC50.lsq/xSmooth)^Hill.lsq)
#	
#	ySmooth.MLE <- 100 / (1 + (IC50.MLE/xSmooth)^Hill.MLE)
#	
#	lines(xSmooth,ySmooth.leastsquares,col = "red")
#	lines(xSmooth,ySmooth.MLE,col = "blue")
#	grid()
#	
#	legendstring1 <- sprintf("least-squares (pIC50 = %0.5g, Hill = %0.5g)",pIC50.lsq,Hill.lsq)
#	legendstring2 <- sprintf("MLE (pIC50 = %0.5g, Hill = %0.5g)",pIC50.MLE,Hill.MLE)
#	
#	legend("bottomright",
#			legend = c(legendstring1,legendstring2),
#			col = c("red","blue"),
#			lty = c("solid","solid"))
#}
#dev.off()

#Run the MCMC chain for a much shorter number of samples, to try and get a good value for the jump parameter.
#====
print("Calculating jump parameter")

#iman#
#mleVector=startingVector
#burnin.jumpcalc <- burnin#10000
#nsamp.jumpcalc <-  nsamp*thin #10000
jumpFactor <- 1.00 #Brad changed to 1.00 for publication purposes.      ----------------------------------
upperAcceptableLimit <- 0.90
lowerAcceptableLimit <- 0.30
acceptanceRatio <- Inf
failureCounter <- 0
failureCounterMax <- 200
#print(sprintf("Burn-In for jump calculation = %d.",burnin.jumpcalc))
#print(sprintf("# of samples for jump calculation = %d.",nsamp.jumpcalc))
#print(sprintf("Lower acceptable limit = %0.5g.",lowerAcceptableLimit))
#print(sprintf("Upper acceptable limit = %0.5g.",upperAcceptableLimit))
#print(sprintf("Starting; failure counter = %d. Maximum allowed failures = %d.",failureCounter,failureCounterMax))
#print("Using initial parameters:")
#print(mleVector)
#print("Gradient at initial point:")
#iman 
#Iman
addburn1 = 20000
j_switch=1
cov_switch=0
Last_step_cova=matrix()
addburn <- burnin
bypass <- FALSE
jump_opt=mleVector * jumpFactor # Iman for initial condition
if(bypass == FALSE)
{
	converged <- FALSE
	while((acceptanceRatio < lowerAcceptableLimit | acceptanceRatio > upperAcceptableLimit) & (failureCounter <= failureCounterMax) | (converged == FALSE))
	{
#		print("j_swich cov_swich")
#		print(j_switch)
#		print(cov_switch)
		
		#print("Input parameters for jump calculation:")
		#print(mleVector * jumpFactor)
		#print(sprintf("Inputting jump factor = %0.5g",jumpFactor))
		print("jumpFactor befor go to MCMC")
		print(jumpFactor)
		
		myChain <- modMCMC(f = compute_likelihood8,
				p = mleVector,
				inputFrame = myframe,
				negation = -2, #iman ?
				printData = FALSE,
				pIC50LocationShift = -4, #iman ?
				priorParameters = priorParameters,
				jump = jump_opt ,#*j_swich+Last_step_cova*cov_swich,
				lower = LB,
				upper = UB,
				wvar0 = 0,
				var0 = NULL,
				updatecov = 10,
				ntrydr = 2,
				burninlength = burnin*j_switch,
				niter = burnin + nsamp*thin*j_switch,
				outputlength = nsamp)
		
		
		#iman 
		
#		pdf(sprintf("visual_test1_%d.pdf",seednum),width=9, height=6)
#		traceplot(myChain)
#		dev.off
		
#		png(sprintf("visual_test%d_%d.png",seednum,burnin))
#		plot(myChain)
#		dev.off()
#		png(sprintf("visual_par_test%d_%d.png",seednum,burnin))
#		plot(myChain$pars)
#		dev.off()
		
		
		nAccepted <- myChain$naccepted
		acceptanceRatio <- nAccepted/(burnin + nsamp * thin*j_switch)
		#acceptanceRatio <- calculate_real_acceptance_rate(mcmc_chain = myChain, chain_length = nsamp)
		print(sprintf("Acceptance ratio = %0.5g",acceptanceRatio))
		print(sprintf("Total aterations after jump part= %d",Total_itr))
		if(acceptanceRatio > upperAcceptableLimit)
		{
			failureCounter <- failureCounter + 1
			print(sprintf("Failures = %d.",failureCounter))
			print("acceptanceRatio too high; increasing jump factor.")
			jumpFactor <- jumpFactor * 1.05
			jump_opt=mleVector * jumpFactor
			print("jumpFactor in 1.05")
			print(jumpFactor)
		}
		
		if(acceptanceRatio < lowerAcceptableLimit)
		{
			failureCounter <- failureCounter + 1
			print(sprintf("Failures = %d.",failureCounter))
			print("acceptanceRatio too low; decreasing jump factor.")
			jumpFactor <- jumpFactor * 0.95
			#jumpFactor <- jumpFactor * 1/2
	        jump_opt=mleVector * jumpFactor
			print("jumpFactor in 0.95")
			print(jumpFactor)
		}
		
		if((acceptanceRatio < upperAcceptableLimit) & (acceptanceRatio > lowerAcceptableLimit))
		{
			print("Successfully found an estimate for the jump parameter.")
			print("Now testing if the Geweke diagnostic is OK.")
			ti=ti+1
			geweke_diagnostic_test <- geweke.diag(as.mcmc(myChain$pars))
			Total_itr=burnin+Total_itr
#			png(sprintf("geweke_visual_test%d_%d.png",seednum,burnin))
#			geweke.plot(as.mcmc(myChain$pars))
#			dev.off()
# 			if (ti>2){
#			png(sprintf("geweke_visual_test%d_%d_%d.png",seednum,burnin,ti))
#			geweke.plot(as.mcmc(myChain$pars))
#			dev.off()
#		    }
			
#			png(sprintf("1geweke_visual_test%d_%d_%d.png",seednum,burnin,ti))
#			geweke.plot(myChain, .1 , .5)
#			dev.off()
#			
#			pdf(sprintf("1geweke_visual_test%d_%d_%d.pdf",seednum,burnin,ti))
#			geweke.plot(myChain, .1 , .5)
#			dev.off()
			
			if(all(abs(geweke_diagnostic_test[["z"]]) < qnorm(0.975)) | Total_itr>Itr_UL)
			{
				converged <- TRUE
				print("Geweke diagnostic: OK!")
				print(sprintf("Using jumpFactor = %0.5g",jumpFactor))
				print(geweke_diagnostic_test[["z"]])
			}
			
#			if(Total_itr>Itr_UL)
#			{
#				converged <- TRUE
#				print("Upper Limit of Iteratoins Achived ")
#				#print(sprintf("Using jumpFactor = %0.5g",jumpFactor))
#			}
			
			
			else
			{
				#print("Geweke diagnostic: NOT OK! Increasing burnin...")
				#burnin <- burnin + addburn
				
				print("Geweke diagnostic: NOT OK! Increasing burnin...")
				#burnin <- burnin + addburn
				burnin <-  addburn1
				#iman
		     #   if (j_swich==1){
				nro=nrow(myChain$pars)
				mleVector=myChain$pars[nro,]
			#}
			#if (j_swich==0){
			#	j_off=0
			#	cov_on=1
				jump_opt=cov(myChain$pars)
				j_switch=0
				#print(jump_opt)
			#}
				
			}
		}
		
		
	}
	
	
	
	if(failureCounter >= failureCounterMax)
	{
		dir.create("logfiles/failed_runs",showWarnings = FALSE)
		write("",sprintf("logfiles/failed_runs/failed_drug_%s",targetDrug))
		stop("Could not find a good initial jump factor.")
	}
}

#if(nrow(myJump) == ncol(myJump))
#{
#	jump_cov_msg <- sprintf("Using covariance matrix for initial jump with %d rows and %d columns.",nrow(myJump),ncol(myJump))
#    print(jump_cov_msg)
#}

#Running test section
#print("Running test section.")
#print("mleVector:")
#print(mleVector)
#print("jumpFactor:")
#print(jumpFactor)
#myChain_test <- modMCMC(f = compute_likelihood8,
#		p = mleVector,
#		inputFrame = myframe,
#		negation = -2,
#		printData = FALSE,
#		pIC50LocationShift = -4,
#		priorParameters = priorParameters,
#		jump = mleVector * jumpFactor,
#		lower = LB,
#		upper = UB,
#		wvar0 = 0,
#		var0 = NULL,
#		updatecov = 10,
#		ntrydr = 2,
#		burninlength = burnin.jumpcalc,
#		niter = burnin.jumpcalc + nsamp.jumpcalc)
#print("passed test section.")


print(warnings())
print("starting MCMC chain")
print(sprintf("Using jumpFactor = %0.5g",jumpFactor))
myJump <-  jumpFactor * mleVector
print("Using jump vector of jumpFactor * mleVector:")
print(myJump)
#ACTUAL CHAIN
#myChain_final <- modMCMC(f = compute_likelihood8,
#		p = mleVector,
#		inputFrame = myframe,
#		negation = -2,
#		printData = FALSE,
#		pIC50LocationShift = -4,
#		priorParameters = priorParameters,
#		jump = mleVector * jumpFactor,
#		lower = LB,
#		upper = UB,
#		wvar0 = 0,
#		var0 = NULL,
#		updatecov = 10,
#		ntrydr = 2,
#		burninlength = burnin,
#		niter = burnin + nsamp*thin)
myChain_final <- myChain

#myChain <- modMCMC(f = compute_likelihood8,
#		p = mleVector,
#		inputFrame = myframe,
#		negation = -2,
#		printData = FALSE,
#		pIC50LocationShift = -4,
#		priorParameters = priorParameters,
#		#jump = rep(jumpFactor,length(mleVector)),
#		jump = mleVector * jumpFactor,
#		lower = LB,
#		upper = UB,
#		wvar0 = 0,
#		var0 = NULL,
#		updatecov = 10,
#		ntrydr = 2,
#		burninlength = burnin.jumpcalc,
#		niter = burnin.jumpcalc + nsamp.jumpcalc)

colnames(myChain_final$pars) <- names(mleVector)
location_tag <- paste0(targetDrug,"_",targetChannel,"_pIC50_location")

pIC50_location_samples <- myChain_final$pars[,location_tag]
IC50_samples <- 10^(6 - pIC50_location_samples)

cmaxData <- read.csv("../../Input/Drug Cmaxes/newCiPA_wendy_noflec_nometh.csv")
cmaxes <- cmaxData$therapeutic
drugnames <- cmaxData$drug
names(cmaxes) <- drugnames
freecmax <- cmaxes[targetDrug] / 1e9
#iman
freecmax=freecmax*10^6

myquantiles <- quantile(IC50_samples,c(.025,.50,.975))

safetyMarginBounds <- myquantiles / freecmax
dir.create("/scratch/bradley.ridder/Hierarchical_Modeling_redone_concs_sites_assumed_confirmed/Output/FOLDER_1/",showWarnings = FALSE)
output.filename.csv <- sprintf("/scratch/bradley.ridder/Hierarchical_Modeling_redone_concs_sites_assumed_confirmed/Output/FOLDER_1/%s_%s.csv",targetDrug,seednum)
print(sprintf("Successful run. Writing %s.",output.filename.csv))
write.csv(data.frame(myChain_final$pars),output.filename.csv,row.names = FALSE)

thinnedFrame <- data.frame(myChain_final$pars)
thinnedFrame <- thinnedFrame[seq(thin,nrow(thinnedFrame),by = thin),]

dir.create("/scratch/bradley.ridder/Hierarchical_Modeling_redone_concs_sites_assumed_confirmed/Output/FOLDER_2/",showWarnings = FALSE)
output.filename.csv.thinned <- sprintf("/scratch/bradley.ridder/Hierarchical_Modeling_redone_concs_sites_assumed_confirmed/Output/FOLDER_2/%s_%s.csv",targetDrug,seednum)
print(sprintf("Successful run. Writing %s.",output.filename.csv.thinned))
write.csv(thinnedFrame,output.filename.csv.thinned,row.names = FALSE)

if(Total_itr>Itr_UL) #check and store the cases that reach the maximum iteration
{
dir.create("/scratch/bradley.ridder/Hierarchical_Modeling_redone_concs_sites_assumed_confirmed/Output/FOLDER_3/",showWarnings = FALSE)
output.filename.csv.thinned <- sprintf("/scratch/bradley.ridder/Hierarchical_Modeling_redone_concs_sites_assumed_confirmed/Output/FOLDER_3/%s_%s.csv",targetDrug,seednum)
print(sprintf("Successful run. Writing %s.",output.filename.csv.thinned))
write.csv(thinnedFrame,output.filename.csv.thinned,row.names = FALSE)
}


print(safetyMarginBounds)