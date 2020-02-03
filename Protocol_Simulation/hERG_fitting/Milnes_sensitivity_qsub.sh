#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=2G
#$ -S /bin/sh
#$ -j y
#$ -pe thread 1
#$ -t 1-28          
#$ -N ms_qsub
#$ -o logfiles_ms_qsub

# Bradley Ridder
# 5/4/2018
# Create Milnes sensitivity plots in parallel.


source /home/kelly.chang/R/source.sh

DRUGNAMES=(dofetilide bepridil sotalol quinidine \
           cisapride terfenadine ondansetron chlorpromazine \
           verapamil ranolazine mexiletine diltiazem \
           ibutilide disopyramide vandetanib azimilide \
           droperidol domperidone clarithromycin risperidone pimozide clozapine \
           metoprolol loratadine tamoxifen nifedipine nitrendipine \
           astemizole)
           
#IDX=$(((SGE_TASK_ID-1)%$WORKERPERDRUG+1))
#SAMPLEPERWORKER=2000/$WORKERPERDRUG
#ISTART=$(((IDX-1)*$SAMPLEPERWORKER+1))
#ISTOP=$((IDX*$SAMPLEPERWORKER))
#ISTR="$ISTART"-"$ISTOP"
#IDX2=$(((SGE_TASK_ID-IDX)/$WORKERPERDRUG))
#DRUG=${DRUGNAMES[IDX2]}
#DOSE="1-4"
#IDX=$((SGE_TASK_ID-1))

DRUG=${DRUGNAMES[$SGE_TASK_ID]}

# Rscript AP_simulation.R -d "$DRUG" -x "$DOSE" -i "$ISTR" -c "data/newCiPA.csv" -b "../hERG_fitting/results/" -m "../IC50Brad/results/" >& logfiles/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID"

Rscript Milnes_sensitivity.R -d "$DRUG" >& logfiles_ms_qsub/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID"