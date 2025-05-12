#!/bin/bash
#SBATCH --job-name=k2_lizard
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

module load kraken2
# kraken2-build --standard --db kraken2_std_db

cd /ibex/user/majnouym/cs249_a2/lizard_asembly/kraken

DB=kraken2_std_db
THREADS=64

for fq in \
  /ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_rev.fastq.gz \
  /ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz
do
  base=$(basename "$fq" .fastq.gz)
  kraken2 \
    --db /ibex/user/majnouym/cs249_a2/lizard_asembly/kraken/kraken2_std_db \
    --threads "$THREADS" \
    --report "${base}.kraken.report.txt" \
    --output "${base}.kraken.out.txt" \
    "$fq"
done