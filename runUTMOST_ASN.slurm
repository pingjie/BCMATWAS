#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=5-12:00:00
#SBATCH --mail-user=jie.ping.1@vumc.org
#SBATCH --mail-type=ALL
#SBATCH --job-name=UTM_ASN
#SBATCH --array=0-212
#SBATCH --output=/nobackup/sbcs/pingj2/EURASN_TWAS/scripts/log/UTMOST_ASN.%A_%a.log
#SBATCH --error=/nobackup/sbcs/pingj2/EURASN_TWAS/scripts/err/UTMOST_ASN.%A_%a.err

nodeName=`hostname`
echo "hostname is: " $nodeName
echo "slurm job ID is: " $SLURM_JOBID
echo "SLURM_ARRAY_JOB_ID is: " $SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

module load PLINK/1.9b_5.2 GCC/10.2.0 OpenMPI/4.0.5 R/4.0.5

race=ASN
while read line
do
  Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/runUTMOST.R ${line} ${race}
done < /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/genes_${race}_utmost/${race}_UTMOST.genes_${SLURM_ARRAY_TASK_ID}
