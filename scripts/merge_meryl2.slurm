#!/bin/bash
#SBATCH --job-name=merge2
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=256G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

module load merqury

cd /ibex/user/majnouym/cs249_a2/lizard_asembly/merqury_2

meryl union-sum output /ibex/user/majnouym/cs249_a2/new_assembly/meryl/lizard.meryl /ibex/user/majnouym/cs249_a2/new_assembly/meryl/lizard_rev.k21.meryl /ibex/user/majnouym/cs249_a2/new_assembly/meryl/lizard_seq.k21.meryl