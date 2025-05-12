#!/bin/bash
#SBATCH --job-name=merq_lizard
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=400G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

cd /ibex/user/majnouym/cs249_a2/lizard_asembly/merqury_2

module load merqury

$MERQURY/merqury.sh /ibex/user/majnouym/cs249_a2/lizard_asembly/merqury_2/lizard.meryl /ibex/user/majnouym/cs249_a2/lizard_asembly/hifiasm/lizard_asm.fasta merq