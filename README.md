# Comprehensive *in vitro* Proarrhythmia Assay (CiPA)<br/> *in silico* model
R code used to perform Lab-specific Calibration and Validation Strategy  for the CiPA *in silico* model


## Authors
Xiaomei Han, Mohammadreza Samieegohar, Zhihua Li
(modified by Aroli Marcellinus for the custom needs of CML)

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


## DISCLAIMER
This research was partially supported by the Ministry of Food and Drug Safety  (24212MFDS266), National Research Foundation of Korea (NRF) under the Basic Science Research Program (2022R1A2C2006326), IITP-Innovative Human Resource Development for Local Intellectualization program grant funded by the Korea government(MSIT)(IITP-2025-RS-2020-II201612), and the Technology development Program(RS-2024-00510494) funded by the Ministry of SMEs and Startups(MSS, Korea).

## References
*Han, X., Samieegohar, M., Ridder, B., J., Wu, W. W,. Randolph, A., Tran, P., Sheng, J., Stoelzle-Feix, S., Brinkwirth, N.,Rotordam, M., G., Becker, N., Friis, S., Rapedius, M., Goetze, A., T., Strassmaier, Tim., Okeyo, G., Kramer, J., Kuryshev, Y., Wu, C., Strauss, D., G.,    Li, Z., et al. "A Lab-specific Calibration and Validation Strategy for Implementing Proarrhythmia Risk Prediction Models: A Case Study of CiPA". (submitted)
