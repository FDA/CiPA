#!/bin/bash

# File:         run_hERG_fit_boot.sh
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Bash script to fit hERG model to bootstrap samples of Milnes
#               protocol data.
#

# change CMA-ES default hyperparameters
POP_SIZE="80"
STOPTOL="0.001"

DRUGNAMES=(dofetilide bepridil sotalol quinidine cisapride terfenadine ondansetron chlorpromazine verapamil ranolazine mexiletine diltiazem)

for DRUG in ${DRUGNAMES[@]}; do
    for IDX in `seq 0 99`; do
        ISTART=$((IDX*20+1))
        ISTOP=$(((IDX+1)*20))
        ISTR="$ISTART"-"$ISTOP"
        Rscript hERG_fitting.R -d "$DRUG" -i "$ISTR" -l "$POP_SIZE" -t "$STOPTOL" > logfiles/"$DRUG"."$ISTR" 2>&1
    done
done
