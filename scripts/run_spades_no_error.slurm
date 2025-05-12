#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16GB
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err

module load spades/3.15.5

cd /ibex/user/majnouym/cs249_a2/

spades.py -s data/synthetic_dataset/reads/no_error_ont_hq_50x.fastq -o outputs/task_1_3_4/ont/no_error --phred-offset 33