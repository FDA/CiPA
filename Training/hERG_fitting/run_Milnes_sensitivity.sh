#!/bin/bash

# File:         run_Milnes_sensitivity.sh
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Bash script to perform sensitivity analysis with bootstrap
#               samples.
#

DRUGNAMES=(dofetilide bepridil sotalol quinidine cisapride terfenadine ondansetron chlorpromazine verapamil ranolazine mexiletine diltiazem)

for DRUG in ${DRUGNAMES[@]}; do
    Rscript Milnes_sensitivity.R -d "$DRUG" >& logfiles/"$DRUG".Milnes
done
