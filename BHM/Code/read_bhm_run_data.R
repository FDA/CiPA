read_bhm_run_data <- function(file_path)
{
  run_parameters <- read.table(file_path,
                               header = TRUE,
                               stringsAsFactors = FALSE,
                               colClasses = c("character","character"))
  
  drug_names <- run_parameters[run_parameters$parameter_type == "drug",
                               "parameter_value"]
  
  seed_values <- run_parameters[run_parameters$parameter_type == "seed",
                                "parameter_value"]
  seed_values <- as.integer(seed_values)
  
  channel_name <- run_parameters[run_parameters$parameter_type == "channel",
                                 "parameter_value"]
  
  total_final_samples <- run_parameters[run_parameters$parameter_type == "total_final_samples",
                                        "parameter_value"]
  total_final_samples <- as.integer(total_final_samples)
  
  burn_length <- run_parameters[run_parameters$parameter_type == "burn_length",
                                "parameter_value"]
  burn_length <- as.integer(burn_length)
  
  thinning_rate <- run_parameters[run_parameters$parameter_type == "thinning_rate",
                                  "parameter_value"]
  thinning_rate <- as.integer(thinning_rate)
  
  
  phase_1_directory <- run_parameters[run_parameters$parameter_type == "phase_1_directory",
                                  "parameter_value"]
  
  phase_2_directory <- run_parameters[run_parameters$parameter_type == "phase_2_directory",
                                      "parameter_value"]
  
  successful_mcmc_chain_output_folder <- run_parameters[run_parameters$parameter_type == "successful_mcmc_chain_output_folder",
                                      "parameter_value"]
  
  exceeded_max_iterations_chain_output_folder <- run_parameters[run_parameters$parameter_type == "exceeded_max_iterations_chain_output_folder",
                                                        "parameter_value"]
  
  top_level_prior_parameters <- run_parameters[run_parameters$parameter_type == "top_level_prior_parameters",
                                                                "parameter_value"]
  
  initial_jump <- run_parameters[run_parameters$parameter_type == "initial_jump",
                                "parameter_value"]
  
  initial_jump <- as.numeric(initial_jump)
  
  output_list <- list()
  output_list$drug_names <- drug_names
  output_list$seed_values <- seed_values 
  output_list$channel_name <- channel_name
  output_list$total_final_samples <- total_final_samples
  output_list$burn_length <- burn_length
  output_list$thinning_rate <- thinning_rate
  output_list$phase_1_directory <- phase_1_directory
  output_list$phase_2_directory <- phase_2_directory
  output_list$successful_mcmc_chain_output_folder <- successful_mcmc_chain_output_folder
  output_list$exceeded_max_iterations_chain_output_folder <- exceeded_max_iterations_chain_output_folder
  output_list$top_level_prior_parameters <- top_level_prior_parameters
  output_list$initial_jump <- initial_jump
  return(output_list)
}