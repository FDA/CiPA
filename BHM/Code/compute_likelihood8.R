source("dllogis.R")
source("dshiftedgamma.R")
#The likelihood function used for the BHM Markov-chain Monte Carlo sampling.
#Takes in dose-response data for an ion channel and a given drug.
#negation is a final multiplier on the likelihood function. For use with the FME package, it will be -2.
#printData is for debugging purposes; leave as FALSE.
#pIC50LocationShift follows what was done in Johnstone et al. 2013. Leave as -4.
#priorParameters are from the Elkins et al. data set used in the Johnstone paper.
compute_likelihood8 <- function(thetaVector,
                                inputFrame,
                                negation = 1,
                                printData = FALSE,
                                pIC50LocationShift = -4,
                                priorParameters)
{
  drugVector <- inputFrame$drug
  concVector <- inputFrame$conc
  channelVector <- inputFrame$channel
  siteVector <- inputFrame$sitename
  blockVector <- inputFrame$block
  
  if(length(unique(drugVector)) != 1)
  {
    stop("drugVector has multiple drugs in it; likelihood function assumes only a single drug.")
  }
  
  if(length(unique(channelVector)) != 1)
  {
    stop("channelVector has multiple channels in it; likelihood function assumes only a single channel.")
  }
  
  logLikelihood <- 0 
  for(i in 1:length(drugVector))
  {
    myDrug <- drugVector[i]
    myChannel <- channelVector[i]
    mySite <- siteVector[i]
    actualBlock <- blockVector[i]
    x <- concVector[i] / 1e9
    x=x*10^6
    pIC50tag <- paste0(myDrug,"_",myChannel,"_",mySite,"_pIC50")
    pIC50 <- thetaVector[pIC50tag]
    
    Hilltag <- paste0(myDrug,"_",myChannel,"_",mySite,"_Hill")
    Hill <- thetaVector[Hilltag]
    
    sigma <- thetaVector["globalsigma"]
    
    IC50 <- 10^(6 - pIC50)
    
    calculatedBlock <- 100/(1 + (IC50/x)^Hill)
    
    logLikelihood <- logLikelihood + dnorm(actualBlock,
                                           mean = calculatedBlock,
                                           sd = sigma,
                                           log = TRUE)
    
    
  }
  
  #Now include other likelihood aspects from the hierarchical model.
  getHill <- grepl("_Hill$",names(thetaVector))
  hillCoefficients <- thetaVector[getHill]
  
  #Get the Hill shape and scale parameters.
  hillShapeName <- paste0(unique(drugVector),"_",unique(channelVector),"_Hill_shape")
  hillScaleName <- paste0(unique(drugVector),"_",unique(channelVector),"_Hill_scale")
  
  for(i in 1:length(hillCoefficients))
  {
    logLikelihood <- logLikelihood + dllogis(xdata = hillCoefficients[i],
                                             shapeParameter = thetaVector[hillShapeName],
                                             scaleParameter = thetaVector[hillScaleName],
                                             log = TRUE)
  }
  
  #Now do pIC50.
  getpIC50 <-  grepl("_pIC50$",names(thetaVector))
  pIC50parameters <- thetaVector[getpIC50]
  
  #Get the pIC50 location and scale parameters.
  pIC50LocationName <- paste0(unique(drugVector),"_",unique(channelVector),"_pIC50_location")
  pIC50ScaleName    <- paste0(unique(drugVector),"_",unique(channelVector),"_pIC50_scale")
  
  for(i in 1:length(pIC50parameters))
  {
    logLikelihood <- logLikelihood + dlogis(x = pIC50parameters[i],
                                            location = thetaVector[pIC50LocationName],
                                            scale = thetaVector[pIC50ScaleName],
                                            log = TRUE)
  }
  
  #And then add in the contributions from the prior distributions for the hyperparameters.
  #Hill Shape
  logLikelihood <- logLikelihood + dgamma(thetaVector[hillShapeName],
                                          shape = priorParameters$Hill.shape.k,
                                          scale = priorParameters$Hill.shape.theta,
                                          log = TRUE)
  #Hill Scale
  logLikelihood <- logLikelihood + dgamma(thetaVector[hillScaleName],
                                          shape = priorParameters$Hill.scale.k,
                                          scale = priorParameters$Hill.scale.theta,
                                          log = TRUE)
  
  #pIC50 location
  logLikelihood <- logLikelihood + dshiftedgamma(thetaVector[pIC50LocationName],
                                                 shapeParameter = priorParameters$pIC50.location.k,
                                                 scaleParameter = priorParameters$pIC50.location.theta,
                                                 shiftParameter = pIC50LocationShift,
                                                 log = TRUE)
  
  #pIC50 scale
  logLikelihood <- logLikelihood + dgamma(thetaVector[pIC50ScaleName],
                                          shape = priorParameters$pIC50.scale.k,
                                          scale = priorParameters$pIC50.scale.theta,
                                          log = TRUE)
  
  #globalsigma scale
  logLikelihood <- logLikelihood + dgamma(thetaVector["globalsigma"],
                                          shape = priorParameters$globalsigma.k,
                                          scale = priorParameters$globalsigma.theta,
                                          log = TRUE)
  
  if(printData == TRUE)
  {
    if(negation < 0)
    {
      print(sprintf("-L = %0.5g.",negation*logLikelihood))
    }
    else
    {
      print(sprintf("L = %0.5g.",negation*logLikelihood))
    }
  }
  return(negation*logLikelihood)
}