#!/bin/bash

# File:         run_AP_uncertainty.sh
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Bash script to run AP simulations for the 12 CiPA training
#               drugs using uncertainty inputs (i.e. drug parameters from
#               sampling distributions).
#

DRUGNAMES=(dofetilide bepridil sotalol quinidine cisapride terfenadine ondansetron chlorpromazine verapamil ranolazine mexiletine diltiazem)

for DRUG in ${DRUGNAMES[@]}; do
    for IDX in `seq 0 99`; do
        ISTART=$((IDX*20+1))
        ISTOP=$(((IDX+1)*20))
        ISTR="$ISTART"-"$ISTOP"
        Rscript AP_simulation.R -d "$DRUG" -i "$ISTR" > logfiles/"$DRUG"."$ISTR" 2>&1
    done
done
