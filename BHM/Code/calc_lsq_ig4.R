source("dllogis.R")
source("parameter_subfit_pIC50.R")
#Calculates a good initial guess from the experimental data.
#This initial guess will be used for maximum-likelihood estimation.
#The initial guess will probably be very close to the MLE result.
calc_lsq_ig4 <- function(masterframe)
{
  #print("Getting combinations.")
  allCombinations <- expand.grid(unique(masterframe$drug),
                                 unique(masterframe$channel),
                                 unique(masterframe$sitename),
                                 stringsAsFactors = FALSE)
  #print("Got combinations.")
  vectorNames <- vector()
  for(i in 1:nrow(allCombinations))
  {
    drugName <- allCombinations[i,1]
    channelName <- allCombinations[i,2]
    site <- allCombinations[i,3]
    subframe <- subset(masterframe,
                       subset = (drug == drugName) & (channel == channelName) & (sitename == site))
    if(nrow(subframe) == 0)
    {
      next
    }
    assembledName <- paste0(drugName,"_",channelName,"_",site)
    
    vectorNames <- c(vectorNames,assembledName)
  }
  vectorNames <- unique(vectorNames)
  vnpIC50 <- paste0(vectorNames,"_pIC50")
  vnHill <- paste0(vectorNames,"_Hill")
  fullVectorNames <- c(vnpIC50,vnHill)
  
  parameterVector <- numeric(length(fullVectorNames))
  names(parameterVector) <- fullVectorNames
  failures <- 0
  failedFits <- data.frame(drug = character(),channel = character(),sitename <- character())
  
  for(i in 1:nrow(allCombinations))
  {
    drugName <- allCombinations[i,1]
    channelName <- allCombinations[i,2]
    site <- allCombinations[i,3]
    subframe <- subset(masterframe,
                       subset = (drug == drugName) & (channel == channelName) & (sitename == site))
    if(nrow(subframe) == 0)
    {
      #print(c(drugName,channelName,site))
      next
    }
    #print(sprintf("Performing fit #%d out of %d",i,nrow(allCombinations)))
    #print(sprintf("drug = %s, channel = %s, site = %s.",drugName,channelName,site))
    outputSolution <- parameter_subfit_pIC50(subframe)
    
    if(outputSolution$optimizationData$convergence != 0)
    {
      failures <- failures + 1
      myFailure <- data.frame(drug = drugName,channel = channelName,sitename = site)
      failedFits <- rbind(failedFits,myFailure)
    }
    
    targetName_pIC50 <- paste0(drugName,"_",channelName,"_",site,"_pIC50")
    parameterVector[targetName_pIC50] <- outputSolution$pIC50
    
    targetName_Hill <- paste0(drugName,"_",channelName,"_",site,"_Hill")
    parameterVector[targetName_Hill] <- outputSolution$Hill
  }
  
  parameterVector <- c(parameterVector,0)
  tempStore <- names(parameterVector)
  tempStore[length(tempStore)] <- "globalsigma"
  names(parameterVector) <- tempStore
  
  allCombinations_hyper <- expand.grid(unique(masterframe$drug),
                                       unique(masterframe$channel),
                                       stringsAsFactors = FALSE)
  
  hyperNames <- vector()
  for(i in 1:nrow(allCombinations_hyper))
  {
    drugName <- allCombinations_hyper[i,1]
    channelName <- allCombinations_hyper[i,2]
    subframe <- subset(masterframe,
                       drug == drugName,
                       channel == channelName)
    if(nrow(subframe) == 0)
    {
      next
    }
    assembledName <- paste0(drugName,"_",channelName)
    #Each drug-channel pIC50 gets its own set of hyperparameters.
    #pIC50 follows a logistic distribution. Needs mu and s parameters (location and scale).
    #Hill follows a log-logistic distribution. Needs alpha and beta parameters (shape and scale).
    #globalsigma follows a gamma distribution. Needs k and theta parameters (shape and scale).
    hyperNames <- c(hyperNames,assembledName)
  }
  
  hyperNames <- unique(hyperNames)
  
  hyperParameters.pIC50.location <- paste0(hyperNames,"_pIC50_location")
  hyperParameters.pIC50.scale <- paste0(hyperNames,"_pIC50_scale")
  
  hyperParameters.Hill.shape <- paste0(hyperNames,"_Hill_shape")
  hyperParameters.Hill.scale <- paste0(hyperNames,"_Hill_scale")
  
  fullHyperNames <- c(hyperParameters.pIC50.location,
                      hyperParameters.pIC50.scale,
                      hyperParameters.Hill.shape,
                      hyperParameters.Hill.scale)
  
  hyperParameterVector <- numeric(length(fullHyperNames))
  names(hyperParameterVector) <- fullHyperNames
  
  completeParameterVector <- c(parameterVector,hyperParameterVector)
  
  #Compute globalsigma estimate.
  mydata <- parameter_subfit_pIC50(masterframe)
  completeParameterVector["globalsigma"] <- mydata$sigma
  
  mlefun.pIC50 <- function(xpars,pIC50data)
  {
    return(-sum(dlogis(pIC50data,location = xpars[1],scale = xpars[2],log = TRUE)))
  }
  
  mlefun.Hill <- function(xpars,Hilldata)
  {
    return(-sum(dllogis(Hilldata,shape = xpars[1],scale = xpars[2],log = TRUE)))
  }
  
  getpIC50 <-  grepl("_pIC50$",names(completeParameterVector))
  pIC50dataInput <- completeParameterVector[getpIC50]
  
  mle.output.pIC50 <- optim(par = c(mean(pIC50dataInput),sd(pIC50dataInput)),
                            fn = mlefun.pIC50,
                            pIC50data = pIC50dataInput,
                            method = "L-BFGS-B",
                            lower = c(-Inf,0.001))
  
  pIC50LocationName <- paste0(unique(masterframe$drug),"_",unique(masterframe$channel),"_pIC50_location")
  pIC50ScaleName    <- paste0(unique(masterframe$drug),"_",unique(masterframe$channel),"_pIC50_scale")
  
  completeParameterVector[pIC50LocationName] <- mle.output.pIC50$par[1]
  completeParameterVector[pIC50ScaleName]    <- mle.output.pIC50$par[2]
  
  
  getHill <- grepl("_Hill$",names(completeParameterVector))
  hilldataInput <- completeParameterVector[getHill]
  hillShapeName <- paste0(unique(masterframe$drug),"_",unique(masterframe$channel),"_Hill_shape")
  hillScaleName <- paste0(unique(masterframe$drug),"_",unique(masterframe$channel),"_Hill_scale")
  
  
  
  mle.output.Hill <- optim(par = c(mean(hilldataInput),sd(hilldataInput)),
                           fn = mlefun.Hill,
                           Hilldata = hilldataInput,
                           method = "L-BFGS-B",
                           lower = c(0.001,0.001))
  
  completeParameterVector[hillShapeName] <- mle.output.Hill$par[1]
  completeParameterVector[hillScaleName]    <- mle.output.Hill$par[2]
  
  return(completeParameterVector)
}