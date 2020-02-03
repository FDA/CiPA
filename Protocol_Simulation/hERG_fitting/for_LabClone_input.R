drugNames <- list.files("results")
for(myDrug in drugNames)
{
	inputString <- sprintf("results/%s/boot_pars.csv",myDrug)
	myFrame <- read.csv(inputString)
	writeName <- sprintf("results/%s_boot_pars.csv",myDrug)
	write.csv(myFrame,writeName,row.names = FALSE)
}