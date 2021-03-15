#!/usr/bin/env bash
#SBATCH --output=test.out
#SBATCH --job-name=solarPCE_parser
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blp62@pitt.edu
#SBATCH --time=0-01:00:00

# Load modules
module purge
module use /ihome/ghutchison/dlf57/software/privatemodules/smp
module load anaconda3

# run script
python solarPCE_parser.py

