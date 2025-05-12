#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

module load merqury

cd /ibex/user/majnouym/cs249_a2/lizard_asembly/merqury_2

meryl union-sum output /ibex/user/majnouym/cs249_a2/lizard_asembly/merqury_2/lizard.meryl /ibex/user/majnouym/cs249_a2/lizard_asembly/merqury_2/lizard_rev.k21.meryl /ibex/user/majnouym/cs249_a2/lizard_asembly/merqury_2/lizard_seq.k21.meryl