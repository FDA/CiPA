#$ -cwd
#$ -l h_rt=10:00:00
#$ -l h_vmem=2G
#$ -S /bin/sh
#$ -j y
#$ -pe thread 1
#$ -t 1-28
#$ -N ss_30sr_ALL
#$ -o logfiles

source /projects/mikem/applications/centos7/R-3.6.0/set-env.sh

DRUGNAMES=(bepridil dofetilide sotalol quinidine \
           cisapride terfenadine ondansetron chlorpromazine \
           verapamil ranolazine mexiletine diltiazem \
           ibutilide disopyramide vandetanib azimilide \
           droperidol domperidone clarithromycin astemizole risperidone pimozide clozapine \
            metoprolol loratadine tamoxifen nitrendipine nifedipine)

IDX=$((SGE_TASK_ID-1))
DRUG=${DRUGNAMES[IDX]}

Rscript Bmax_Rscript.R \
-d "$DRUG" \
-f data/drug_block_hERG_30_second_ramp_protocol.csv \
-c "hERG" \
>& logfiles/"$DRUG"