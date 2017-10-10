#!/bin/bash

# File:         run_IC50_mcmc.sh
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Bash script to run Markov-chain Monte Carlo simulations for
#               the 12 CiPA training drugs.
#

DRUGNAMES=(dofetilide bepridil sotalol quinidine cisapride terfenadine ondansetron chlorpromazine verapamil ranolazine mexiletine diltiazem)

mkdir logfiles
for DRUG in ${DRUGNAMES[@]}; do
    Rscript IC50_mcmc.R -d "$DRUG" >& logfiles/"$DRUG"
done
