#!/bin/bash

# File:         run_herg_bootstrap.sh
# Author:       Kelly Chang
# Date:         Sep 2017
# Modified by:	Aroli Marcellinus
# Date:		May 2025
# Version:      1.0
# 
# Description:  Bash script to demonstrate the use of hERG bootstrap fitting
#               code. This code has been modified for the usage of hERG bootstrap
#

mkdir results
mkdir figs

if [ "$#" -lt 3 ]; then
    echo "Error: Provide experimental time-series file, drug name, and number of samples that want to be bootstrapped and the sample size"
    exit 1
fi


#--- create bootstrap samples
HERG_FILE=$1
DRUG=$2
SAMPLE_SIZE=$3

# Create 2000 bootstrap samples. The output of this command is used later for
# fitting.

Rscript generate_bootstrap_samples.R -d "$DRUG" -f "$HERG_FILE" -n "$SAMPLE_SIZE"

#--- compile model
cd models
R CMD SHLIB hergmod.c
cd ..


#--- fitting specifications

# Number of cores to use for parallel evaluation (optional, should not affect
# reproducibility)
NCORES=$(nproc --all)

#--- cmaes hyperparameters (optional, omit for defaults)
# Set population size larger to reduce the chance of getting a local minimum
# and for quicker convergence.
POP_SIZE="80"

# Set stopping tolerance higher than default (error changes little below
# this tolerance).
STOPTOL="0.001"

# Set maximum number of generations low for this example only. The default
# number of generations is recommended--omit this option (-m) to use the
# default.
MAX_GEN="10"

#--- fit optimal parameters (output required for bootstrap fitting)
Rscript hERG_fitting.R -f "$HERG_FILE" -d "$DRUG" -c "$NCORES" -l "$POP_SIZE" -m "$MAX_GEN" -t "$STOPTOL"

#--- fit 10 bootstraps
# Quick example with only 10 bootstraps (more should be performed for the
# real analysis).
BOOTNUM="1-$SAMPLE_SIZE"
Rscript hERG_fitting.R -f "$HERG_FILE" -d "$DRUG" -i "$BOOTNUM" -c "$NCORES" -l "$POP_SIZE" -m "$MAX_GEN" -t "$STOPTOL"

#--- combine results, plot sampling distributions, and calculate CI
Rscript process_boot_results.R -d "$DRUG"

#--- get sensitivity of block during Milnes protocol
Rscript Milnes_sensitivity.R -d "$DRUG" -f "$HERG_FILE"
