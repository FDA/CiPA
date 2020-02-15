#!/bin/sh
#$ -cwd
#$ -pe thread 1
#$ -j y
#$ -N hERG_boot
#$ -l s_rt=24:00:00
#$ -R y
#$ -l h_rt=24:00:00
#$ -t 1-28000          #num. drugs(28) x WORKERPERDRUG(1000)
#$ -o logfiles

DRUGNAMES=(dofetilid ebepridil sotalol quinidine cisapride terfenadine ondansetron chlorpromazine verapamil ranolazine mexiletine diltiazem ibutilide disopyramide vandetanib azimilide droperidol domperidone clarithromycin risperidone pimozide clozapine metoprolol loratadine tamoxifen nifedipine nitrendipine astemizole)


NDRUG=${#DRUGNAMES[@]}
IDX=$(((SGE_TASK_ID-1)/NDRUG))
BOOTNUM=$((IDX+1))
IDX2=$((SGE_TASK_ID-1-IDX*NDRUG))
DRUG=${DRUGNAMES[IDX2]}

POP_SIZE=80
STOPTOL=0.001
MAXIT=100

Rscript hERG_fitting.R -d "$DRUG" -i "$BOOTNUM" -l "$POP_SIZE" -m "$MAXIT" -t "$STOPTOL" >& logfiles/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID"


