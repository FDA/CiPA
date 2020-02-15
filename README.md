# Comprehensive *in vitro* Proarrhythmia Assay (CiPA)<br/> *in silico* model
R code used to perform Lab-specific Calibration and Validation Strategy  for the CiPA *in silico* model

## Authors
Xiaomei Han, Mohammadreza Samieegohar, Zhihua Li

## Requirements
This code was developed with R version 3.5 and uses the following packages:
* optparse (version 1.4.4)
* deSolve (version 1.14)
* cmaes (version 1.0-11)
* FME (version 1.3.5)
* rms (4.5-0)
* ggplot2 (version 1.14)

## Uncertainty characterization
[hERG_fitting/](hERG_fitting/) contains code to perform uncertainty characterization for drug binding kinetics of the human Ether-Ã -go-go-Related Gene (hERG) channel gating model ([README.md](hERG_fitting/README.md)).
Note: compared to the previous/original hERG data in CiPAORdv1.0 (Li, Zhihua, et al.* 2019)( https://github.com/FDA/CiPA/tree/Model-Validation-2018/hERG_fitting/data) a new filtering step was applied to select only cells with less than 20% background current. More information is provided in the supplementary materials of the paper “A Lab-specific Calibration and Validation Strategy for Implementing Proarrhythmia Risk Prediction Models: A Case Study of CiPA” (submitted).

[chantest_Hill_fitting/](chantest_Hill_fitting/) contains code to perform uncertainty characterization of dose-response curves for three ionic currents (ICaL, INaL and INa). 28 CiPA drugs data provides by Charles Rivers Laboratories. ([README.md]( chantest_Hill_fitting/README.md)).
Note: compared to the previous/original dose response data from Charles Rivers Laboratories in CiPAORdv1.0 (Li, Zhihua, et al.* 2019) (https://github.com/FDA/CiPA/tree/Model-Validation-2018/Hill_Fitting/data), a new filtering step was applied to replace cells with less than 60% block at highest concentration. This filtering step has resulted in the replacement of the original verapamil data on ICaL by new data from the same company. More information is provided in the supplementary materials of the paper “A Lab-specific Calibration and Validation Strategy for Implementing Proarrhythmia Risk Prediction Models: A Case Study of CiPA” (submitted).

[nanion_Hill_fitting/](nanion_Hill_fitting/) contains code to perform uncertainty characterization of dose-response curves for three ionic currents (ICaL, INaL and INa). 27 CiPA drugs data provides by Nanion Technologies GmbH.([README.md](nanion_Hill_fitting/README.md)).

Note: In both nanion and chantest Hill_fittings, the Markov Chain Monte Carlo (MCMC) algorithm was updated from the previous version (https://github.com/FDA/CiPA/blob/Model-Validation-2018/Hill_Fitting/IC50_mcmc.R). The updated code is under the Hill fitting folders for Nanion and Chantest. More information is provided in the supplementary material of the paper “A Lab-specific Calibration and Validation Strategy for Implementing Proarrhythmia Risk Prediction Models: A Case Study of CiPA” (submitted).


## Uncertainty propagation
[chantest_AP_simulation/](chantest_AP_simulation/) contains code to propagate uncertainty in drug effects to action potential (AP) simulations. Results from [hERG_fitting/](hERG_fitting/) and [chantest_Hill_fitting/](chantest_Hill_fitting/) are used as model inputs ([README.md](chantest_AP_simulation/README.md)).

[nanion8rand_AP_simulation/](nanion8rand_AP_simulation/) contains code to propagate uncertainty in drug effects to action potential (AP) simulations. Results from [hERG_fitting/](hERG_fitting/) and [nanion_Hill_fitting/](nanion_Hill_fitting/) are used as model inputs ([README.md](nanion_Hill_fitting/README.md)).


[sens/]( sens/) contains codes to perform sensitivity analysis of all 28 drugs where a small perturbation was applied to each drug’s torsade metric score, and the corresponding change in the classification thresholds are recorded. This step filtered all 28 drugs to a list of top influential drugs that will be selected as the candidate list to screen for calibration drugs. Results from [chantest_AP_simulation/](chantest_AP_simulation/) are used as inputs ([README.md]([sens/](sens/))).

[LOOCV/](LOOCV/) contains codes to perform leave-one-out-cross-validation and compute the 8 performance metric values using calibration drugs..

## DISCLAIMER
This software and documentation were developed by the authors in their capacities as Oak Ridge Institute for Science and Education (ORISE) research fellows at the U.S. Food and Drug Administration (FDA).

FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.  Further, FDA makes no representations that the use of the Software will not infringe any patent or proprietary rights of third parties.   The use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions.

## References
*Han, X., Samieegohar, M., Ridder, B., J., Wu, W. W,. Randolph, A., Tran, P., Li, Z., et al. "A Lab-specific Calibration and Validation Strategy for Implementing Proarrhythmia Risk Prediction Models: A Case Study of CiPA". (submitted)
