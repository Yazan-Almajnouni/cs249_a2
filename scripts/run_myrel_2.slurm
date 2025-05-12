#!/bin/bash
#SBATCH --job-name=meryl_lizard
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=256G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL


DIR=$1
OUT=$2
READ=$3


cd $DIR

module load merqury

echo "Running meryl count k=21 output ${OUT} ${READ}"
meryl count k=21 threads=128 output "${OUT}" "${READ}"