#Source the files here.
source("is_between.R")
source("is_csv.R")
source("all_field_names_match.R")
source("load_hERG_fit_parameters.R")
source("get_start_time.R")
source("run_simulation.R")
source("markov_states_sum_to_one.R")
source("markov_states_greater_than_zero.R")
source("generate_events_voltage_dataframe.R")
source("find_maximum_open_probability.R")
source("blockade_calculation.R")
source("drug_blockade_experiment.R")

#load libraries.
library(deSolve)
library(optparse)

#Specify command-line arguments.
#--- specify command line arguments
#-d,--drug
#-c,--concfile
#-m,--modelname
#-e,--modelpath
#-p,--parampath
#-s,--statepath
#-h,--hERGpath
#-T,--temperature
#-o,--outputdir
#-V,--voltagefunction
#-P,--voltagefunctionparameters

for(my_drug in drug_list)
{
  for(target_index in 1:num_pars)
  {
    
  }
}


args<-parse_args(parser)

drugNames <- args$drug
modelname <- args$modelname
pars <- read.table(args$parampath)
parNames <- pars$V1
pars <- pars$V2
names(pars) <- parNames
pars["T"] <- args$temperature
myLoopConcsFrame <- read.csv(args$concfile,stringsAsFactors = FALSE)
dyn.load(args$modelpath)
initstates <- read.table(args$statepath)
stateNames <- initstates$V1
initstates <- initstates$V2
names(initstates) <- stateNames
voltagefunction <- args$voltagefunction
voltagefunctionparameters <- args$voltagefunctionparameters
source(paste0(voltagefunction,".R",collapse=""))
voltagefunction <- match.fun(voltagefunction)
voltage_program_data <- voltagefunction(read.csv(voltagefunctionparameters,
                                            stringsAsFactors = FALSE))
voltage_program <- voltage_program_data$voltage_program
timepoints <- voltage_program_data$timepoints
print(sprintf("Length of voltage program = %0.5g milliseconds.",max(timepoints)))
print(sprintf("Maximum value of voltage program = %0.5g mV.",max(voltage_program)))
print(sprintf("Minimum value of voltage program = %0.5g mV.",min(voltage_program)))
print(sprintf("Run start time: %s.",date()))
print(sprintf("Running drug %s.",args$drug))
print(sprintf("Using concentration file %s.",args$concfile))
print(sprintf("Using modelname = %s.",args$modelname))
print(sprintf("Using model path = %s.",args$modelpath))
print(sprintf("Using parameter path %s.",args$parampath))
print(sprintf("Using initial state path %s.",args$statepath))
print(sprintf("Using hERG parameter path %s.",args$hERGpath))
print(sprintf("Running at temperature %d Celsius.",args$temperature))
print(sprintf("Using hERG sample %d.",args$targetindex))
print(sprintf("Results .csv files will write to %s.",args$outputdir))
print(sprintf("Using voltage function: %s.",args$voltagefunction))
print(sprintf("Using voltage function parameters: %s.",voltagefunctionparameters))
print("Concentrations used (nM):")
print(myLoopConcsFrame)

dir.create(args$outputdir,recursive = TRUE)

for(myDrug in drugNames)
{
	myhERGpath <- sprintf("%s/%s_boot_pars.csv",args$hERGpath,myDrug)
	print(sprintf("Loading hERG parameter data from %s.",myhERGpath))
	hERGparams <- read.csv(myhERGpath)
	myOutputFrame <- data.frame(drug = character(),
			conc = numeric(),
			units = character(),
			block = numeric(),
			channel = character())
	
	for(i in 1:nrow(hERGparams))
	{
		pars["Kmax"] <- hERGparams$Kmax[targetindex]
		pars["Ku"] <- hERGparams$Ku[targetindex]
		pars["n"] <- hERGparams$n[targetindex]
		pars["halfmax"] <- hERGparams$halfmax[targetindex]
		pars["Vhalf"] <- hERGparams$Vhalf[targetindex]
		
		for(myConc in myLoopConcsFrame[[myDrug]])
		{
			print(sprintf("Running hERG sample %d for drug %s at concentration %0.5g nM.",targetindex,myDrug,myConc))
			
			myBlockData <- drug_blockade_experiment(modelname = modelname,
					parameterVector = pars,
					initialStates = initstates,
					drugConc = myConc,
					voltage_program = voltage_program,
					timepoints = timepoints,
					breakTolerance = 0.01,
					maxPulses = 100L,
					minimumPulses = 10L,
					mydrug = myDrug,
					myTemperature = args$temperature)
			
			print(sprintf("Calculated block = %0.5g percent.",100*myBlockData$fractionBlockade))
			placeholderFrame <- data.frame(drug = myDrug,
					conc = myConc,
					units = "nM",
					block = myBlockData$fractionBlockade,
					channel = "hERG")
			
			myOutputFrame <- rbind(myOutputFrame,placeholderFrame)
		}
	}
	dir.create(sprintf("%s/%s",args$outputdir,myDrug))
	outputFilename <- sprintf("%s/%s/%s_%d.csv",args$outputdir,myDrug,myDrug,targetindex)
	print(sprintf("Writing output file %s.",outputFilename))
	write.csv(myOutputFrame,outputFilename,row.names = FALSE)
}
print(sprintf("Completed running drug %s.",myDrug))
print(sprintf("Finish time: %s.",date()))