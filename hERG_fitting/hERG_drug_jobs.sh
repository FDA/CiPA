#!/bin/sh
#$ -cwd
#$ -pe thread 2
#$ -j y
#$ -N hERG_drug
#$ -l s_rt=24:00:00
#$ -R y
#$ -l h_rt=24:00:00
#$ -t 1
#$ -o logfiles


#DRUGNAMES=(dofetilide bepridil sotalol quinidine cisapride terfenadine ondansetron chlorpromazine verapamil ranolazine mexiletine diltiazem)
DRUGNAMES=(ibutilide2)
DRUGNAMES=(cisapride24)

IDX=$((SGE_TASK_ID-1))
DRUG=${DRUGNAMES[IDX]}
POP_SIZE=80
STOPTOL=0.001
MAXIT=4000

Rscript hERG_fitting.R -d "$DRUG" -c 2 -l "$POP_SIZE" -m "$MAXIT" -t "$STOPTOL" >& logfiles/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID"
