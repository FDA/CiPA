#!/bin/bash

DRUGNAMES=(dofetilide bepridil sotalol quinidine cisapride terfenadine ondansetron chlorpromazine verapamil ranolazine mexiletine diltiazem)

for DRUG in ${DRUGNAMES[@]}; do
    Rscript Milnes_sensitivity.R -d "$DRUG" >& logfiles/"$DRUG".Milnes
done
