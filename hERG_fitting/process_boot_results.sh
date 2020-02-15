#!/bin/sh
#$ -cwd
#$ -pe thread 1
#$ -j y
#$ -N hERG_boot
#$ -l s_rt=24:00:00
#$ -R y
#$ -l h_rt=24:00:00
#$ -t 1-15
#$ -o logfiles

source /home/kelly.chang/R/source.sh

#DRUGNAMES=(dofetilide bepridil sotalol quinidine cisapride terfenadine ondansetron chlorpromazine verapamil ranolazine mexiletine diltiazem)
DRUGNAMES=(astemizole azimilide clarithromycin clozapine disopyramide domperidone droperidol ibutilide loratadine metoprolol nifedipine pimozide risperidone tamoxifen vandetanib)
IDX=$((SGE_TASK_ID-1))
DRUG=${DRUGNAMES[IDX]}

Rscript process_boot_results.R -d "$DRUG" >& logfiles/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID"
