#!/bin/bash
#SBATCH --job-name=ABM_Test
#SBATCH --partition=savio
#SBATCH --account=fc_mgenes
#SBATCH --qos=savio_normal
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=dougherty.eric@gmail.com
module load r
R CMD BATCH ABM_SAVIO_Test.R
rm -rf .RData