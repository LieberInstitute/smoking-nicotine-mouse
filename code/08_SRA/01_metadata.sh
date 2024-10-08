#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=3G
#SBATCH --job-name=01_metadata
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ../../processed-data/08_SRA/01_metadata.txt
#SBATCH -e ../../processed-data/08_SRA/01_metadata.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.4

## List current modules for reproducibility
module list

Rscript 01_metadata.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/