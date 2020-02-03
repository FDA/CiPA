#Loads the dose-response data from .csv files into a data frame that can be used for evaluating the MCMC likelihood
#function.
#
#trim is the fraction of the data to load; used for debugging purposes, but otherwise should be 1.0.
#
#siteExclusionList is a vector of strings describing the sites to exclude from loading.
#This works fine as-is, or the same effect can be achieved by moving data not to be used to
#a folder other than ../../Input/phase1_four_concentrations and
#../../Input/phase2_four_concentrations/
#
#allowedChannels is a vector of strings, which can be hERG, INa, INaL, or ICaL.
#
#targetDrug is the drug for which data is being loaded across all sites. In this work, it is
#the 28 CiPA drugs.
makeMasterFrameFaster <- function(trim = 1.0,
                                  phase_1_directory = "../../Input/phase1_four_concentrations/",
                                  phase_2_directory = "../../Input/phase2_four_concentrations/",
                                  siteExclusionList = "none",
                                  allowedChannels,
                                  targetDrug)
{
  phase1.dir <- phase_1_directory
  phase2.dir <- phase_2_directory
  
  #Print inputs.
  #print(sprintf("Trim = %0.5f",trim))
  #print("Excluded sites:")
  #print(siteExclusionList)
  #print("Allowed channels:")
  #print(allowedChannels)
  #print(sprintf("Target drug = %s.",targetDrug))
  
  phaseList <- c(phase1.dir,phase2.dir)
  
  channelFolderNames <- c("cav","herg","late_nav","peak_nav")
  names(channelFolderNames) <- c("ICaL","hERG","INaL","INa")
  channels <- channelFolderNames[allowedChannels]
  
  masterframe <- data.frame(drug = character(),
                            conc = numeric(),
                            units = character(),
                            block = numeric(),
                            channel = character(),
                            sitename = character())
  
  drugExclusionList <- c("nifedipine_1_phase2",
                         "nifedipine_2_phase2",
                         "verapamil_2_phase2",
                         "verapamil_1_phase2",
                         "nifedipinemetrioncontrol",
                         "cisapride_1_phase2",
                         "cisapride_2_phase2",
                         "dofetilidemetrioncontrol",
                         "terfenadinemetrioncontrol",
                         "verapramilmetrioncontrol",
                         "ranolazinemetrioncontrol",
                         "amitriptylinemetrioncontrol")
  print("Excluded drugs:")
  print(drugExclusionList)
  for(k in 1:length(phaseList))
  {
    
    for(i in 1:length(channels))
    {
      filePath <- paste0(phaseList[k],channels[i])
      siteFolderInfo <- list.files(filePath)
      nFiles <- length(siteFolderInfo)
      for(j in 1:nFiles)
      {
        getsitename <- strsplit(siteFolderInfo[j],"_")
        sitename <- getsitename[[1]][1]
        if(sitename %in% siteExclusionList)
        {
          next
        }
        specificPath <- paste0(filePath,"/",siteFolderInfo[j])
        rawdata <- read.csv(specificPath)
        print(specificPath)
        rawdata$sitename <- rep(sitename,nrow(rawdata))
        print(names(rawdata))
        print(names(masterframe))
        masterframe <- rbind(masterframe,rawdata)
      }
    }
  }
  isBadDrug <- (masterframe$drug %in% drugExclusionList) | 
    !(masterframe$channel %in% allowedChannels) |
    (masterframe$drug != targetDrug)
  masterframe <- masterframe[!isBadDrug,]
  
  #Trim the masterframe to compute faster gradient.
  masterframe <- masterframe[1:floor(trim*nrow(masterframe)),]
  #Unique sites loaded.
  print("Unique sites loaded:")
  print(sprintf("%s",unique(masterframe$sitename)))
  return(masterframe)
}