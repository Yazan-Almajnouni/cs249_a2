#!/bin/bash
#SBATCH --job-name=busco_lizard
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

cd /ibex/user/majnouym/cs249_a2

module load busco

busco \
  -i lizard_asembly/hifiasm/lizard_asm.fasta \
  -o lizard_busco \
  --out_path lizard_asembly/busco \
  --auto-lineage \
  -m genome \
  -c 32