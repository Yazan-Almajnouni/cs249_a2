#!/bin/bash
#SBATCH --job-name=merq_lizard
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=400G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

cd /ibex/user/majnouym/cs249_a2/new_assembly/meryl

module load merqury

$MERQURY/merqury.sh /ibex/user/majnouym/cs249_a2/new_assembly/meryl/lizard.meryl /ibex/user/majnouym/cs249_a2/new_assembly/filtered_reads/hifiasm/new_asm.fasta merq