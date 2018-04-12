#!/bin/bash

# File:         run_compute_TdP_error.sh
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Bash script to perform logistic regression on AP simulation
#               results.
#

# classification with fixed-input simulation results
Rscript compute_TdP_error.R > logfiles/compute_fixed_error 2>&1

# classification with uncertainty-input simulation results
Rscript compute_TdP_error.R -u > logfiles/compute_uncertainty_error 2>&1

# classification with uncertainty-input results with omitted current effects
Rscript compute_TdP_error.R -u -o > logfiles/compute_uncertainty_dropcurrent_error 2>&1
