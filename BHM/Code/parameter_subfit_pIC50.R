#Used for calculation of initial guess.
#HillLower and HillUpper should be left as 0.5 and 2.0.
parameter_subfit_pIC50 <- function(inputFrame,
                                   HillLower = 0.5,
                                   HillUpper = 2.0)
{
  if(any(inputFrame$units != "nM"))
  {
    stop("Concentration units of inputFrame are not in nM.")
    print(inputFrame)
  }
  
  print(sprintf("Lower permissible Hill coefficient limit = %0.5g.",HillLower))
  print(sprintf("Upper permissible Hill coefficient limit = %0.5g.",HillUpper))
  
  #Convert nanomolar to micromolar.
  xData <- inputFrame$conc / 1e9
  xData=xData*10^6
  yData <- inputFrame$block

  lsq_fun <- function(thetaParameters,xData,yData)
  {
    IC50 <- 10^(6 - thetaParameters[1])
    Hill <- thetaParameters[2]
    HillEquation <- 100/(1 + (IC50/xData)^Hill)
    SSE <- sum((yData - HillEquation)^2)
    return(SSE)
  }
  
  #Initial guesses.
  #Hill guess of 1 is reasonable.
  #IC50 can be guessed from the data itself.
  xPos <- which.min((yData - 50)^2)
  pIC50IG <- 6 - log10(xData[xPos])
  
  #iterate over pIC50 initial guess of +/- 3 decades until convergence is achieved.
  igvector.lefthalf.pIC50 <- seq(pIC50IG - 3,pIC50IG,length.out = 50)
  igvector.righthalf.pIC50 <- seq(pIC50IG,pIC50IG + 3,length.out = 50)
  igvector <- c(igvector.lefthalf.pIC50,
                igvector.righthalf.pIC50)
  
  bestFunctionValue <- Inf
  bestSolution <- NULL
  for(i in 1:length(igvector))
  {
    myIG <- igvector[i]
    mySolution <- optim(c(myIG,1),
                        lsq_fun,
                        xData = xData,
                        yData = yData,
                        method = "L-BFGS-B",
                        lower = c(-Inf,HillLower),
                        upper = c(Inf,HillUpper))

    if((mySolution$value < bestFunctionValue) & (mySolution$convergence == 0))
    {
      bestFunctionValue <- mySolution$value
      bestSolution <- mySolution
    }
  }
  outputSolution <- list()
  
  if(is.null(bestSolution) == TRUE)
  {
    print("Failed to converge when computing initial guess.")
    print(sprintf("drug = %s.",inputFrame$drug[1]))
    print(sprintf("channel = %s.",inputFrame$channel[1]))
    print(sprintf("sitename = %s.",inputFrame$sitename[1]))
    outputSolution$pIC50 <- pIC50IG
    outputSolution$Hill <- 1
    
    failedIC50 <- 10^(6 - pIC50IG)
    failedHill <- 1
    yFail <- 100/(1 + (failedIC50/xData)^failedHill)
    sigma <- sd(yFail - yData)
    outputSolution$sigma <- sigma
    outputSolution$optimizationData <- mySolution
    return(outputSolution)
  }
  mySolution <- bestSolution

  IC50best <- 10^(6 - mySolution$par[1])
  Hillbest <- mySolution$par[2]
  outputSolution$pIC50 <- mySolution$par[1]
  outputSolution$Hill <- mySolution$par[2]
  
  #Compute sigma of errors.
  yMean <- 100/(1 + (IC50best/xData)^Hillbest)
  residuals <- yMean - yData
  sigma <- sd(residuals)
  
  outputSolution$sigma <- sigma
  outputSolution$optimizationData <- mySolution
  
  return(outputSolution)
}