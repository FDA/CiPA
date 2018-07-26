#$ -cwd
#$ -l h_rt=10:00:00
#$ -l h_vmem=2G
#$ -S /bin/sh
#$ -j y
#$ -pe thread 1
#$ -t 1-28
#$ -N IC50_mcmc
#$ -o logfiles

#Change above "-t 1-28" to either "-t 1-12" for TRAINING SET or  "-t 1-16" for VALIDATION SET.

#Uncomment either the TRAINING SET or the VALIDATION SET, but not both.

#TRAINING SET
#DRUGNAMES=(dofetilide bepridil sotalol quinidine cisapride terfenadine ondansetron chlorpromazine verapamil ranolazine mexiletine diltiazem)

#VALIDATION SET
DRUGNAMES=(droperidol domperidone disopyramide clarithromycin ibutilide metoprolol azimilide risperidone pimozide loratadine tamoxifen clozapine nifedipine nitrendipine vandetanib)

IDX=$((SGE_TASK_ID-1))
DRUG=${DRUGNAMES[IDX]}

Rscript IC50_mcmc.R -d "$DRUG" -f "../../Input/HTS_DrugBlock/HTS_validationdrug_block.csv" >& logfiles/"$DRUG"
