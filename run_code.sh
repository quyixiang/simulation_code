#!/bin/bash
#SBATCH --job-name=simulation
#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 03-00:00:00
#SBATCH --output=./Report/slurm-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=yqu@unc.edu

# cd is really important
cd /proj/cosd_lab/dwuWrite/yixiang/simulation_code

source /nas/longleaf/home/yixiang/.bashrc

module load r
# Rscript real_data_code/Yixiang_create_marginal.R
Rscript simulation_Poisson.R