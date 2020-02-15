#$ -cwd
#$ -l h_rt=10:00:00
#$ -l h_vmem=2G
#$ -S /bin/sh
#$ -j y
#$ -pe thread 1
#$ -t 1-28
#$ -N IC50_mcmc
#$ -o logfiles

#source /home/kelly.chang/R/source.sh

DRUGNAMES=(\dofetilide bepridil sotalol quinidine \
           cisapride terfenadine ondansetron chlorpromazine \
           verapamil ranolazine mexiletine diltiazem \
           ibutilide disopyramide vandetanib azimilide \
           droperidol domperidone clarithromycin astemizole risperidone pimozide clozapine \
            metoprolol loratadine tamoxifen nitrendipine nifedipine\
)

IDX=$((SGE_TASK_ID-1))
DRUG=${DRUGNAMES[IDX]}

Rscript newIC50_mcmc.R -f data/drug_block.csv -d "$DRUG" >& logfiles/"$DRUG"
