read_labclone_run_data <- function(labclone_run_data_list)
{
  run_parameters <- read.table(file_path,
                               header = TRUE,
                               stringsAsFactors = FALSE,
                               colClasses = c("character","character"))
  
  drug_names <- run_parameters[run_parameters$parameter_type == "drug",
                               "parameter_value"]
  
  concfile <- run_parameters[run_parameters$parameter_type == "concfile",
                               "parameter_value"]
  
  concfile <- run_parameters[run_parameters$parameter_type == "concfile",
                             "parameter_value"]
  modelname <- run_parameters[run_parameters$parameter_type == "modelname",
                             "parameter_value"]
  modelpath <- run_parameters[run_parameters$parameter_type == "modelpath",
                             "parameter_value"]
  parampath <- run_parameters[run_parameters$parameter_type == "parampath",
                             "parameter_value"]
  statepath <- run_parameters[run_parameters$parameter_type == "statepath",
                             "parameter_value"]
  hERGpath <- run_parameters[run_parameters$parameter_type == "hERGpath",
                             "parameter_value"]
  temperature <- run_parameters[run_parameters$parameter_type == "temperature",
                             "parameter_value"]
  outputdir <- run_parameters[run_parameters$parameter_type == "outputdir",
                             "parameter_value"]
  voltagefunction <- run_parameters[run_parameters$parameter_type == "voltagefunction",
                             "parameter_value"]
  voltagefunctionparameters <- run_parameters[run_parameters$parameter_type == "voltagefunctionparameters",
                             "parameter_value"]

  
  output_list <- list()
  output_list$drug_names <- drug_names
  output_list$concfile <- concfile
  output_list$modelname <- modelname
  output_list$modelpath <- modelpath
  output_list$parampath <- parampath
  output_list$statepath <- statepath
  output_list$hERGpath <- hERGpath
  output_list$temperature <- temperature
  output_list$outputdir <- outputdir
  output_list$voltagefunction <- voltagefunction
  output_list$voltagefunctionparameters <- voltagefunctionparameters
  return(output_list)
}