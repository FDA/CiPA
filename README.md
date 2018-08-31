# Comprehensive *in vitro* Proarrhythmia Assay (CiPA)<br/> *in silico* model
R code used to perform uncertainty quantification for the CiPA *in silico* model and validate the model with the validation data set.

## Authors
Kelly Chang, Zhihua Li, Bradley J. Ridder

## Requirements
This code was developed with R version 3.3 and uses the following packages:
* optparse (version 1.4.4)
* deSolve (version 1.14)
* cmaes (version 1.0-11)
* FME (version 1.3.5)
* rms (4.5-0)
* ggplot2 (version 1.14)

## Uncertainty characterization
[hERG_fitting/](hERG_fitting/) contains code to perform uncertainty characterization for drug binding kinetics of the human Ether-à-go-go-Related Gene (hERG) channel gating model ([README.md](hERG_fitting/README.md)).

[Hill_fitting/](Hill_fitting/) contains code to perform uncertainty characterization of dose-response curves for six ionic currents (ICaL, INaL, INa, Ito, IKs, and IK1) ([README.md](Hill_fitting/README.md)).

IC50_mcmc_jobs.sh was designed to do training and validation data separately (by providing different input data files).

## Uncertainty propagation
[AP_simulation/](AP_simulation/) contains code to propagate uncertainty in drug effects to action potential (AP) simulations. Results from [hERG_fitting/](hERG_fitting/) and [Hill_fitting/](Hill_fitting/) are used as model inputs ([README.md](AP_simulation/README.md)).

## Work Specific to Validation Paper
Use the hERG_fitting code in the "Training" folder to do hERG fitting on all 28 drugs. Hill_fitting contains four different data sets: 
* manual training
* manual validation
* high-throughput training
* high-throughput validation

Adjust the code in Hill_fitting "IC50_mcmc_jobs.sh" to point at the correct files. An example is given in that code file.

Steps:

* In [/Hill_Fitting](/Hill_Fitting), execute (IC50_mcmc_jobs.sh). Make sure there is a [logfiles/](logfiles/) folder in that folder. A [results/](results/) and [figs/](figs/) folder will appears with results and figures.
* In [/AP_simulation](/AP_simulation), execute (AP_uncertainty_jobs.sh). Adjust AP_uncertainty_jobs.sh to target the results folder you made in the first step.
* In the [/AP_simulation](/AP_simulation) folder, use "Rscript combine_results.R -n 2000" to combine results and produce metrics.rds.

## Paper Citations

The "training" work was reported in the paper:

Li, Z., Dutta, S., Sheng, J., Tran, P.N., Wu, W., Chang, K., Mdluli, T., Strauss, D.G. and Colatsky, T., 2017. Improving the in silico assessment of proarrhythmia risk by combining hERG (human ether-à-go-go-related gene) channel–drug binding kinetics and multichannel pharmacology. Circulation: Arrhythmia and Electrophysiology, 10(2), p.e004628.

The "validation" work was reported in the paper:

Li, Z. , Ridder, B. J., Han, X. , Wu, W. W., Sheng, J. , Tran, P. N., Wu, M. , Randolph, A. , Johnstone, R. H., Mirams, G. R., Kuryshev, Y. , Kramer, J. , Wu, C. , Crumb, W. J. and Strauss, D. G. (2018), Assessment of an In Silico Mechanistic Model for Proarrhythmia Risk Prediction Under the CiPA Initiative. Clin. Pharmacol. Ther. doi:10.1002/cpt.1184

## DISCLAIMER
This software and documentation were developed by the authors in their capacities as Oak Ridge Institute for Science and Education (ORISE) research fellows at the U.S. Food and Drug Administration (FDA).

FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.  Further, FDA makes no representations that the use of the Software will not infringe any patent or proprietary rights of third parties.   The use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions.
