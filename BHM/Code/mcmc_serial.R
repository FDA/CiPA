source("compute_likelihood8.R")
source("calc_lsq_ig4.R")
source("makeMasterFrameFaster.R")
source("read_bhm_run_data.R")
library(FME)

bhm_input_parameters <- read_bhm_run_data(file_path = "../Input/bhm_run_data/bhm_run_data.dat")

seed_values <- bhm_input_parameters$seed_values
drug_names <- bhm_input_parameters$drug_names
seed_values <- bhm_input_parameters$seed_values
channel_name <- bhm_input_parameters$channel_name
total_final_samples <- bhm_input_parameters$total_final_samples 
burn_length <- bhm_input_parameters$burn_length 
thinning_rate <- bhm_input_parameters$thinning_rate 
phase_1_directory <- bhm_input_parameters$phase_1_directory 
phase_2_directory <- bhm_input_parameters$phase_2_directory
successful_mcmc_chain_output_folder <- bhm_input_parameters$successful_mcmc_chain_output_folder
exceeded_max_iterations_chain_output_folder <- bhm_input_parameters$exceeded_max_iterations_chain_output_folder
top_level_prior_parameters <- bhm_input_parameters$top_level_prior_parameters
initial_jump <- bhm_input_parameters$initial_jump

dir.create(successful_mcmc_chain_output_folder,
           showWarnings = FALSE,
           recursive = TRUE)

dir.create(exceeded_max_iterations_chain_output_folder,
           showWarnings = FALSE,
           recursive = TRUE)

for(SEED in seed_values)
{
  for(DRUG in drug_names)
  {
    ti <- 1
    seednum <- SEED
    set.seed(seednum)
    nsamp <- total_final_samples
    burnin <- burn_length
    thin <- thinning_rate
    targetDrug <- DRUG
    targetChannel <- channel_name
    
    print(sprintf("Random seed = %d.",seednum))
    print(sprintf("Channel = %s.",targetChannel))
    print(sprintf("Drug = %s.",targetDrug))
    print(sprintf("nsamp = %d.",nsamp))
    print(sprintf("burnin = %s.",burnin))
    print(sprintf("thin = %s.",thin))
    
    print("Loading dose-response data.")
    myframe <- makeMasterFrameFaster(trim = 1.0,
                                     phase_1_directory = phase_1_directory,
                                     phase_2_directory = phase_2_directory,
                                     siteExclusionList = NULL,
                                     targetChannel,
                                     targetDrug)
    
    #Initial guess calculation by least squares.
    print("Calculating initial guess.")
    startingVector <- calc_lsq_ig4(myframe) 
    print("Success in calculating initial guess.")
    
    myControl <- list()
    myControl$factr <- 1e9
    myControl$maxit <- 1000
    Itr_UL <- 500000
    Total_itr <- 0
    
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
    LB["globalsigma"] <-  0.005
    UB["globalsigma"] <- 50.0
    
    pIC50_scale_tag <- grepl("_pIC50_scale$",g)
    LB[pIC50_scale_tag] <- 0.005
    UB[pIC50_scale_tag] <- Inf
    
    pIC50_location_tag <- grepl("_pIC50_location$",g)
    LB[pIC50_location_tag] <- -3.99
    UB[pIC50_location_tag] <- Inf
    
    Hill_shape_tag <- grepl("_Hill_shape$",g)
    LB[Hill_shape_tag] <- 0.001
    UB[Hill_shape_tag] <- Inf
    
    Hill_scale_tag <- grepl("_Hill_scale$",g)
    LB[Hill_scale_tag] <- 0.001
    UB[Hill_scale_tag] <- Inf
    
    #Load the top-level hierarchy gamma parameters.
    priorParameters <- read.csv(top_level_prior_parameters,
                                stringsAsFactors = FALSE) 
    
    mySolution <- optim(par = startingVector,
                        fn = compute_likelihood8,
                        inputFrame = myframe,
                        negation = -2,
                        printData = FALSE,
                        pIC50LocationShift = -4,
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
    print(sprintf("Length of MLE vector: %d.",length(mleVector)))
    print("startingVector")
    print(startingVector)
    
    print("mleVector")
    print(mleVector)
    
    #Run the MCMC chain for a much shorter number of samples,
    #to try and get a good value for the jump parameter.
    print("Calculating jump parameter")
    jumpFactor <- initial_jump
    upperAcceptableLimit <- 0.90
    lowerAcceptableLimit <- 0.30
    acceptanceRatio <- Inf
    failureCounter <- 0
    failureCounterMax <- 200
    addburn1 <- burnin
    j_switch <- 1
    cov_switch <- 0
    Last_step_cova <- matrix()
    addburn <- burnin
    bypass <- FALSE
    jump_opt <- mleVector*jumpFactor
    
    if(bypass == FALSE)
    {
      converged <- FALSE
      while((acceptanceRatio < lowerAcceptableLimit | acceptanceRatio > upperAcceptableLimit) & (failureCounter <= failureCounterMax) | (converged == FALSE))
      {
        print("jumpFactor calculation before heading to MCMC.")
        print(jumpFactor)
        print(dim(jump_opt))
        myChain <- modMCMC(f = compute_likelihood8,
                           p = mleVector,
                           inputFrame = myframe,
                           negation = -2,
                           printData = FALSE,
                           pIC50LocationShift = -4,
                           priorParameters = priorParameters,
                           jump = jump_opt,
                           lower = LB,
                           upper = UB,
                           wvar0 = 0,
                           var0 = NULL,
                           updatecov = 10,
                           ntrydr = 2,
                           burninlength = burnin*j_switch,
                           niter = burnin + nsamp*thin*j_switch,
                           outputlength = nsamp)
        
        
        
        nAccepted <- myChain$naccepted
        acceptanceRatio <- nAccepted/(burnin + nsamp * thin*j_switch)
        print(sprintf("Acceptance ratio: %0.5g",acceptanceRatio))
        print(sprintf("Total iterations after jump part: %d",Total_itr))
        if(acceptanceRatio > upperAcceptableLimit)
        {
          failureCounter <- failureCounter + 1
          print(sprintf("Failures = %d.",failureCounter))
          print("acceptanceRatio too high; increasing jump factor.")
          jumpFactor <- jumpFactor * 1.05
          jump_opt <- mleVector * jumpFactor
          print(jumpFactor)
        }
        
        if(acceptanceRatio < lowerAcceptableLimit)
        {
          failureCounter <- failureCounter + 1
          print(sprintf("Failures = %d.",failureCounter))
          print("acceptanceRatio too low; decreasing jump factor.")
          jumpFactor <- jumpFactor * 0.95
          jump_opt <- mleVector * jumpFactor
          print(jumpFactor)
        }
        
        if((acceptanceRatio < upperAcceptableLimit) & (acceptanceRatio > lowerAcceptableLimit))
        {
          print("Successfully found an estimate for the jump parameter.")
          print("Now testing if the Geweke diagnostic is OK.")
          ti <- ti + 1
          geweke_diagnostic_test <- geweke.diag(as.mcmc(myChain$pars))
          Total_itr <- burnin+Total_itr
          
          if(all(abs(geweke_diagnostic_test[["z"]]) < qnorm(0.975)) | Total_itr > Itr_UL)
          {
            converged <- TRUE
            print("Geweke diagnostic: OK!")
            print(sprintf("Using jumpFactor = %0.5g",jumpFactor))
            print(geweke_diagnostic_test[["z"]])
          }
          else
          {
            print("Geweke diagnostic: NOT OK! Increasing burnin...")
            burnin <-  addburn1
            nro <- nrow(myChain$pars)
            mleVector <- myChain$pars[nro,]
            jump_opt <- cov(myChain$pars)
            j_switch <- 0
          }
        }
      }
      if(failureCounter >= failureCounterMax)
      {
        dir.create("logfiles/failed_runs",showWarnings = FALSE,recursive = TRUE)
        write("",sprintf("logfiles/failed_runs/failed_drug_%s",targetDrug))
        stop("Could not find a good initial jump factor.")
      }
    }
    myChain_final <- myChain
    myJump <-  jumpFactor * mleVector
    colnames(myChain_final$pars) <- names(mleVector)
    location_tag <- paste0(targetDrug,"_",targetChannel,"_pIC50_location")
    pIC50_location_samples <- myChain_final$pars[,location_tag]
    IC50_samples <- 10^(6 - pIC50_location_samples)
    
    output.filename.csv <- sprintf(paste0(successful_mcmc_chain_output_folder,
                                          "/mcmc_bhm_success_%s_%s.csv"),
                                          DRUG,
                                          SEED)
    
    print(sprintf("Successful run. Writing %s.",output.filename.csv))
    write.csv(data.frame(myChain_final$pars),output.filename.csv,row.names = FALSE)
    
    #Store the cases that reach the maximum iterations
    if(Total_itr > Itr_UL)
    {
      over_iter_file <- sprintf(paste0(exceeded_max_iterations_chain_output_folder,
                                       "/mcmc_bhm_exceeded_max_iter_%d_%s_%s.csv"),
                                       Itr_UL,
                                       DRUG,
                                       SEED)
      
      print(sprintf("Successful run. Writing %s.",over_iter_file))
      write.csv(data.frame(myChain_final$pars),over_iter_file,row.names = FALSE)
    }
    print(sprintf("Complete drug %s for seed %d.",DRUG,SEED))
  }
}