#!/bin/bash

# File:         run_AP_fixed.sh
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Bash script to run AP simulations for the 12 CiPA training
#               drugs using fixed inputs (i.e. optimal drug parameters).
#

DRUGNAMES=(dofetilide bepridil sotalol quinidine cisapride terfenadine ondansetron chlorpromazine verapamil ranolazine mexiletine diltiazem control)

for DRUG in ${DRUGNAMES[@]}; do
    Rscript AP_simulation.R -d "$DRUG" > logfiles/"$DRUG" 2>&1
done
