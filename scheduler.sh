#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=14-00:00:00
#SBATCH --mem=2gb
#SBATCH --job-name=snakemake_bin_wkflw
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#wd = /home/projects-wrighton/NIH_Salmonella/KaiMetaG_20200327/assemble_5.88gbp
eval "$(conda shell.bash hook)"

conda activate metagenomics

#output folder
mkdir -p logs_slurm
#snakemake command
snakemake -s Snakefile --cluster-config cluster_config.yaml --cluster 'sbatch -time={config.time} --mem={config.mem} -o {config.output} -e {config.output} --mail-type {config.email_type} --mail-user {config.email}' -j 2