#!/bin/bash

# File:         run_combine_results.sh
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Bash script to combine AP simulation results.
#

# combine fixed-input simulation results
Rscript combine_results.R

# combine uncertainty-input simulation results
Rscript combine_results.R -n 2000

# combine uncertainty-input results with omitted current effects
Rscript combine_results.R -n 2000 -o
