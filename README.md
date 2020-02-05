# Bayesian Hierarchical Model Code
This code performs Markov-chain Monte Carlo to sample from the posterior distribution suggested by the data supplied by the different laboratory sites.

## BHM Required Packages
This code uses the following R packages:
```
FME
```

To run the code:

```
source("BHM/Code/mcmc_serial.R")
```

## How Parameters are Specified for the BHM Chains

These values are specified in bhm_run_data.dat:

```
(/BHM/)(/Input/bhm_run_data/bhm_run_data.dat)
```

```
parameter_type parameter_value
drug <drugname1>
drug <drugname2>
...
drug <drugnameN>
seed <seedinteger1>
seed <seedinteger2>
...
seed <seedintegerM>
channel hERG
total_final_samples 20000
burn_length 10000
thinning_rate 10
phase_1_directory ../Input/phase1_four_concentrations/
phase_2_directory ../Input/phase2_four_concentrations/
successful_mcmc_chain_output_folder ../Output/successful_mcmc_chain_run_output
exceeded_max_iterations_chain_output_folder ../Output/exceeded_max_iterations_mcmc_chain_run_output
top_level_prior_parameters ../Input/Gamma_Top_Level/Gamma_Top_Level.csv
initial_jump 1.00
```
drug is a drug name
seed is an integer that specifies a random seed to use
channel hERG should not be changed
total_final_samples is the total length of the final chain, but the actual output file will only contain X samples, where X = total_final_samples/thinning_rate.
burn_length is the length of the initial burn-in phase that is discarded. Burn in is used to ignore the transient phase and (theoretically) make sure the chain is only sampled from the stationary distribution.
thinning_rate keep every X'th sample and discards the rest, where X is some integer.
phase_1_directory should not be changed
phase_2_directory should not be changed
successful_mcmc_chain_output_folder is where the output CSV files are written
exceeded_max_iterations_chain_output_folder is where chains that exceeded their iteration allotment are written
top_level_prior_parameters is parameters for the highest level of the hierarchical model. It should not be changed.
initial_jump is how widely the proposal distribution can sample. It is difficult to know ahead of time what is a "good" jump value.
To compensate for this, the code is intended to start at 1.00, and marches down with new initial_jumps until the acceptance ratio is from 0.30 <= r <= 0.90.

The defaults can be used as-is to verify the code runs properly.

# Protocol Simulation Code
The protocol simulation code is a modification and extension of our previous work in <put a reference here>. This code, for a given drug, cycles through all 2,000 parameter sets for the hERG Markov model discussed in <that reference>. 
## Protocol_Simulation Required Packages
This code uses the following R packages:
```
testthat
deSolve
cmaes
optparse
deSolve
cmaes
FME
rms
ggplot2

```

The IKr-dynamic ORd model C code is provided in [models/](models/) and must be compiled:
```
cd BHM_PS_2020
R CMD SHLIB /Protocol_Simulation/Protocol_Simulation/hERG_Model_C_Code/hergmod.c
```

## Protocol_Simulation Unit Tests and testthat
Protocol_Simulation comes with a battery of unit tests for all of the component functions. Running the unit tests successfully is good proof everything is working as it should and repeatable results will be obtained.
To run all of the unit tests, install testthat:

```
install.packages("testthat")
library(testthat)
```

Then, while in the directory (Protocol_Simulation/)(/Protocol_Simulation/)(/Test_Script/), run the command:
```
source("run_labclone_unit_tests.R")
```

This will run all of the ~250 unit tests. All greens means everything is working properly.
Also, the unit tests are the easiest and fastest way to understand how the code works, and how it ought to work.
Those interested in expanding upon this code will likely find them very helpful.

## Running Protocol_Simulation simulations
If all of the unit tests pass, you can run the Protocol_Simulation simulations. To run them, use:

```
cd /Protocol_Simulation/Protocol_Simulation/Run_Script
source("read_labclone_run_data.R")
source("labclone_serial.R")
```

The simulation parameters are specified in:
```
(/Protocol_Simulation/)(/hERG_fitting/)(/results/)
(/Protocol_Simulation/)(/(/Protocol_Simulation/)/)(/Input/)(/hERG_parms/)

```

Results and figures are automatically saved to [Protocol_Simulation/](/Protocol_Simulation/)(/Output/)(/<particular voltage protocol>)

The resulting CSV files will have the drug concentration in one column, and the block in another. This is repeated for all 2000 hERG parameters. For each drug, 2000 CSV files will be created.

## How Parameters are Specified for the Protocol_Simulation Simulations
The values specified in labclone_run_data.dat are:
```
drug <drugname1>
drug <drugname2>
...
drug <drugnameN>
concfile (../Input/)(/input_concentrations/)(/input_concentrations.csv/)
modelname hergmod
modelpath (..Code/)(/hergmod.dll/)
parampath (../Input/)(/hERG_parms/)(/hergmod_pars.txt/)
statepath (../Input/)(/hERG_parms/)(/hergmod_states.txt/)
hERGpath (../Input/)(/hERG_fitting_results_rejection_strategy/)
temperature 37
outputdir (../Output/)(/5_second_ramp_protocol/)
voltagefunction voltage_step_program_list_mode
voltagefunctionparameters (../voltage_program_data_files/)(/5_second_ramp_protocol.csv/)
```

concfile is a CSV containing the concentrations to simulate each drug at.
modelname is always hergmod
modelpath is the compiled .dll files from hergmod.c. It must be compiled on your own computer or cluster before it will run.
parampath is static parameters that should not be changed.
statepath is a .txt file showing the initial states for simulation. This probably should be left as-is.
hERGpath is the path to the individual drug CSV files containing the 2000 hERG parameters.
temperature is in Celsius.
outputdir is where the output files will be written for the Protocol_Simulation simulations.
voltagefunction is how the voltage protocol is specified. The most general choice is voltage_trueAP_program_list_mode, which uses an input CSV files of t and V (time and voltage) to specify any arbitrary voltage protocol.
voltagefunctionparameters is a CSV file that specifies the voltage protocol.

The defaults can be used as-is to verify the code runs properly.
