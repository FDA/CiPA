# Uncertainty characterization of dose-response curves
This code performs uncertainty characterization of dose-response curves (the Hill equation) for ionic current block using a Bayesian inference approach.

## Data format
Patch clamp data should be stored in CSV format with the following headers:

* **drug**: drug name
* **conc**: drug concentration in nM
* **channel**: name of ionic current tested
* **block**: amount of block (%)

To process data used in Li *et al.* 2017 to the correct format, run:

```
cd data
Rscript process_patchclamp_data.R
```

This will create the file [data/drug_block.csv](data/drug_block.csv), containing experimental results for the 12 CiPA training drugs as well as some additional compounds (see [README.md](data/README.md) for details).

## Running the code
This code uses the R packages optparse (version 1.4.4) and FME (version 1.3.5).

To fit bepridil:

```
Rscript IC50_mcmc.R -d bepridil
```

The code attempts to fit an IC50 value and Hill coefficient for each drug-channel pair in the data and then obtains a joint sampling distribution of the parameters using Markov-chain Monte Carlo simulation (MCMC). For data that cannot be fitted, these values are omitted from the output.

New data can be fitted by specifying the data file path and the drug name:

```
Rscript IC50_mcmc.R -d drug1 -f my_data_file.csv
```

By default, 2000 samples are saved from the MCMC run. A different number of samples can be saved:

```
Rscript IC50_mcmc.R -d drug1 -f my_data_file.csv -n 3000
```

For additional help:

```
Rscript IC50_mcmc.R -h
```

## References
* Li, Z., Dutta, S., Sheng, J., Tran, P.N., Wu, W., Chang, K., et al. (2017). Improving the In Silico Assessment of Proarrhythmia Risk by Combining hERG (Human Ether-à-go-go-Related Gene) ChannelDrug Binding Kinetics and Multichannel Pharmacology. Circulation: Arrhythmia and Electrophysiology 10(2), e004628. doi: 10.1161/circep.116.004628.
