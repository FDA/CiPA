# Propagation of uncertainty in drug effects
This code performs uncertainty propagation for the CiPA *in silico* model. Model inputs are the pharmacological parameters for the hERG/IKr Markov model ([README.md](../hERG_fitting/README.md)) and Hill equation parameters for drug block of other ionic currents ([README.md](../Hill_fitting/README.md)). These inputs are used to simulate action potentials (APs) with the optimized IKr-dynamic ORd model (Dutta *et al.* 2017). The primary model output considered here is the qNet metric for proarrhythmia risk described by Dutta *et al.* 2017.

## Running the code
This code uses the following R packages: optparse (version 1.4.4), deSolve (version 1.14), ggplot2 (version 2.2.0), and rms (version 4.5-0).

The IKr-dynamic ORd model C code is provided in [models/](models/) and must be compiled:

```
cd models
R CMD SHLIB newordherg_qNet.c
```

To run simulations, either fixed-input (optimal-fit) parameters or uncertainty-input parameters must be supplied to the model. By default, drug-hERG parameters are located [here](../hERG_fitting/results/) and Hill equation parameters are located [here](../Hill_fitting/results/).

Results and figures are automatically saved to [results/](/results/) and [figs/](figs/), respectively.

[This bash script](run_example.sh) provides a short example of how to run the simulations. The full process is explained below.

### Fixed-input simulations
The control (no drug) AP simulation can be run by calling the script:

```
Rscript AP_simulation.R
```

AP simulations with drug can be run by specifying the drug name (case sensitive) and dose (interpreted as multiples of the therapeutic concentration (nM) listed in [data/CiPA_training_drugs.csv](data/CiPA_training_drugs.csv)):

```
Rscript AP_simulation.R -d dofetilide -x "1-10,15,20,25"
```

New drugs can be simulated if therapeutic (Cmax) values and input parameters are provided at the specified locations:

```
Rscript AP_simulation.R -d drug1 -x "1-10,15,20,25" --cmaxfile="my_cmax_table.csv" --hergpath="path/to/herg/results/" --hillpath="path/to/hill/results/"
```

If the therapeutic concentration is not available in the CSV file specified by "--cmaxfile", doses are interpreted as nM concentrations.

### Uncertainty-input simulations
Samples from uncertainty-input probability distributions can be simulated by specifying sample indices:

```
Rscript AP_simulation.R -d dofetilide -x "1-10,15,20,25" -i "1-2000"
```

The appropriate parameters will be looked up in the default directories ([here](../hERG_fitting/results/) and [here](../Hill_fitting/results/)).

Note that running many AP simulations is computationally intensive, and it is recommended to run them in parallel on a high-performance computing resource. (See [this script](run_AP_uncertainty.sh) for an example of how to split up the simulations.

### Postprocessing
Results for all fixed-input simulations can be combined with:

```
Rscript combine_results.R
```

Results for all uncertainty-input simulations can be combined by specifying the number of samples that were run:

```
Rscript combine_results.R -n 2000
```

The scripts [compute_qNet_CI.R](compute_qNet_CI.R) and [compute_TdP_error.R](compute_TdP_error.R) perform analysis on the combined results and generate figures. With these scripts, TdP risk categories for each drug (0, 1, or 2) are read from a CSV file specified by the "--tdpfile" option ([data/CiPA_training_drugs.csv](data/CiPA_training_drugs.csv) by default):

```
Rscript compute_qNet_CI.R
Rscript compute_TdP_error.R
Rscript compute_TdP_error.R --uncertainty
```

## References
* Dutta, S., Chang, K.C., Beattie, K.A., Sheng, J., Tran, P.N., Wu, W.W., et al. (2017). Optimization of an In silico Cardiac Cell Model for Proarrhythmia Risk Assessment. Frontiers in Physiology 8(616). doi: 10.3389/fphys.2017.00616.
* O'Hara, T., Virag, L., Varro, A., and Rudy, Y. (2011). Simulation of the undiseased human cardiac ventricular action potential: model formulation and experimental validation. PLoS Comput Biol 7(5), e1002061. doi: 10.1371/journal.pcbi.1002061.
