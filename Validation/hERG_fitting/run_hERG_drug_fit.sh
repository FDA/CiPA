#!/bin/bash

# File:         run_hERG_fit_drug.sh
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Bash script to fit hERG model to Milnes protocol data.
#

# change CMA-ES default hyperparameters
POP_SIZE="80"
STOPTOL="0.001"

DRUGNAMES=(dofetilide bepridil sotalol quinidine cisapride terfenadine ondansetron chlorpromazine verapamil ranolazine mexiletine diltiazem)

for DRUG in ${DRUGNAMES[@]}; do
    Rscript hERG_fitting.R -d "$DRUG" -l "$POP_SIZE" -t "$STOPTOL" > logfiles/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID" 2>&1
done
