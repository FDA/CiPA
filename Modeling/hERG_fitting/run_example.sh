#!/bin/bash

# File:         run_example.sh
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Bash script to demonstrate the use of hERG bootstrap fitting
#               code. This code is provided for example only.
#

mkdir logfiles

#--- create bootstrap samples
DRUG="bepridil"

# Create 2000 bootstrap samples. The output of this command is used later for
# fitting.

Rscript generate_bootstrap_samples.R -d "$DRUG"

#--- compile model
cd models
R CMD SHLIB hergmod.c
cd ..

#--- fitting specifications

# Number of cores to use for parallel evaluation (optional, should not affect
# reproducibility)
NCORES="2"

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
Rscript hERG_fitting.R -d "$DRUG" -c "$NCORES" -l "$POP_SIZE" -m "$MAX_GEN" -t "$STOPTOL" > logfiles/"$DRUG" 2>&1

#--- fit 10 bootstraps
# Quick example with only 10 bootstraps (more should be performed for the
# real analysis).
BOOTNUM="1-10"
Rscript hERG_fitting.R -d "$DRUG" -i "$BOOTNUM" -c "$NCORES" -l "$POP_SIZE" -m "$MAX_GEN" -t "$STOPTOL" > logfiles/"$DRUG"."$BOOTNUM" 2>&1

#--- combine results, plot sampling distributions, and calculate CI
Rscript process_boot_results.R -d "$DRUG"

#--- get sensitivity of block during Milnes protocol
Rscript Milnes_sensitivity.R -d "$DRUG" > logfiles/"$DRUG".Milnes 2>&1
