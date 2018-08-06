#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=2G
#$ -S /bin/sh
#$ -j y
#$ -pe thread 1
WORKERPERDRUG=100
#$ -t 1-2800           #num. drugs x WORKERPERDRUG
#$ -N uncertainty
#$ -o logfiles

DRUGNAMES=(\
           dofetilide bepridil sotalol quinidine \
           cisapride terfenadine ondansetron chlorpromazine \
           verapamil ranolazine mexiletine diltiazem \
           ibutilide disopyramide vandetanib azimilide \
           droperidol domperidone clarithromycin risperidone pimozide clozapine \
           metoprolol loratadine tamoxifen nifedipine nitrendipine \
           astemizole\
           )

IDX=$(((SGE_TASK_ID-1)%$WORKERPERDRUG+1))
SAMPLEPERWORKER=2000/$WORKERPERDRUG
ISTART=$(((IDX-1)*$SAMPLEPERWORKER+1))
ISTOP=$((IDX*$SAMPLEPERWORKER))
ISTR="$ISTART"-"$ISTOP"
IDX2=$(((SGE_TASK_ID-IDX)/$WORKERPERDRUG))
DRUG=${DRUGNAMES[IDX2]}

DOSE="1-4"

Rscript AP_simulation.R -d "$DRUG" -x "$DOSE" -i "$ISTR" -c "../data/newCiPA.csv" -m "../Hill_Fitting/results" >& logfiles/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID"
