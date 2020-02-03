#!/bin/bash

DRUGNAMES=(astemizole azimilide bepridil chlorpromazine cisapride clarithromycin clozapine diltiazem disopyramide dofetilide domperidone droperidol ibutilide loratadine metoprolol mexiletine nifedipine nitrendipine ondansetron pimomzide quinidine ranolazine risperidone sotalol tamoxifen terfenadine vandetanib verapamil)

for DRUG in ${DRUGNAMES[@]}; do
    Rscript Milnes_sensitivity.R -d "$DRUG" >& logfiles/"$DRUG".Milnes
done
