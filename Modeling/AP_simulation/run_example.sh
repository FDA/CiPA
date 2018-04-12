#!/bin/bash

# File:         run_example.sh
# Author:       Kelly Chang
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Bash script to demonstrate the use of AP simulation scripts.
#               This code is provided for example only.
#

# run control simulation
Rscript AP_simulation.R >& logfiles/control

# simulate dofetilide with fixed parameters
DRUG="dofetilide"
DOSE="1-10,15,20,25"

Rscript AP_simulation.R -d "$DRUG" -x "$DOSE" >& logfiles/"$DRUG"

# simulate dofetilide with parameters from uncertainty sampling distributions
ISAMP="1-10"
Rscript AP_simulation.R -d "$DRUG" -i "$ISAMP" -x "$DOSE" >& logfiles/"$DRUG"."$ISAMP"

# simulate dofetilide uncertainty sample with individual currents' effects omitted, one at a time
ISAMP="1"
DOSE="1-4"

Rscript AP_simulation.R -d "$DRUG" -o -i "$ISAMP" -x "$DOSE" >& logfiles/"$DRUG"."$ISAMP".drop_currents

# alternatively, specify which current's effects to omit
DROPCURRENT="hERG"
#DROPCURRENT="ICaL"
#DROPCURRENT="INaL"
#DROPCURRENT="INa"
#DROPCURRENT="IKs"
#DROPCURRENT="Ito"
#DROPCURRENT="IK1"

Rscript AP_simulation.R -d "$DRUG" -r "$DROPCURRENT" -i "$ISAMP" -x "$DOSE" >& logfiles/"$DRUG"."$ISAMP".drop_"$DROPCURRENT"
