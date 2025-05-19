#$ -cwd
#$ -l h_rt=10:00:00
#$ -l h_vmem=2G
#$ -S /bin/sh
#$ -j y
#$ -pe thread 1
#$ -t 1-28
#$ -N IC50_mcmc
#$ -o logfiles

if [ "$#" -lt 1 ]; then
    echo "Error: Provide drug name that want to be bootstrapped"
    exit 1
fi

DRUG=$1

Rscript newIC50_mcmc.R -f data/drug_block.csv -d "$DRUG"
