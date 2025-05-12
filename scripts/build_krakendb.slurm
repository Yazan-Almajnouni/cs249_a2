#!/bin/bash
#SBATCH --job-name=buildk
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

module load kraken2

kraken2-build --standard --threads 64 --db /ibex/user/majnouym/cs249_a2/lizard_asembly/kraken/kraken2_std_db 