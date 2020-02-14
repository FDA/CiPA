#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=2G
#$ -S /bin/sh
#$ -j y
#$ -pe thread 1
WORKERPERDRUG=1000
#$ -t 1-27000          #num. drugs x WORKERPERDRUG
#$ -N uncertainty
#$ -o logfiles

#source /home/kelly.chang/R/source.sh
#disopyramidenoICaL
DRUGNAMES=(\
           dofetilide bepridil sotalol quinidine \
           cisapride terfenadine ondansetron chlorpromazine \
           verapamil ranolazine mexiletine diltiazem \
           ibutilide vandetanib azimilide \
           droperidol domperidone clarithromycin risperidone pimozide clozapine \ 
           metoprolol loratadine tamoxifen nifedipine nitrendipine astemizole \ 
           )

IDX=$(((SGE_TASK_ID-1)%$WORKERPERDRUG+1))
SAMPLEPERWORKER=2000/$WORKERPERDRUG
ISTART=$(((IDX-1)*$SAMPLEPERWORKER+1))
ISTOP=$((IDX*$SAMPLEPERWORKER))
ISTR="$ISTART"-"$ISTOP"
IDX2=$(((SGE_TASK_ID-IDX)/$WORKERPERDRUG))
DRUG=${DRUGNAMES[IDX2]}

DOSE="1-4"

Rscript AP_simulation.R -d "$DRUG" -x "$DOSE" -i "$ISTR" -c "data/newCiPA.csv" -m "../nanion_Hill_fitting/results" >& logfiles/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID"
