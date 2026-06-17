#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=2G
#$ -S /bin/sh
#$ -j y
#$ -pe thread 1
#$ -t 1-6
#$ -N fixed
#$ -o logfiles

#source /home/kelly.chang/R/source.sh

DRUGNAMES=(\
           dofetilide bepridil sotalol quinidine \
           cisapride terfenadine ondansetron chlorpromazine \
           verapamil ranolazine mexiletine diltiazem \
           ibutilide disopyramide vandetanib azimilide \
           droperidol domperidone clarithromycin risperidone pimozide clozapine \
           metoprolol loratadine tamoxifen nifedipine nitrendipine \
           astemizole
           )
IDX=$((SGE_TASK_ID-1))
DRUG=${DRUGNAMES[IDX]}

DOSE="0.16,0.32,0.64,0.96,1.28,2.56"

Rscript AP_simulation.R -d "$DRUG" -c "../pureManual_AP_simulation/data/newCiPA.csv" -b "/home/xiaomei.han/cardio/UQ/hERG_fitting/results/" -m "/scratch/xiaomei.han/UQ/chantest_Hill_fitting/results" -x "$DOSE" >& logfiles/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID"
