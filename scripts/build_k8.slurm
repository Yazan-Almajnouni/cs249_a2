#!/bin/bash
#SBATCH --job-name=bldk8
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

module load kraken2

export KRAKEN2_DB=/ibex/user/majnouym/cs249_a2/lizard_asembly/kraken_8

echo "[$(date)] Creating db directory: $KRAKEN2_DB"
mkdir -p $KRAKEN2_DB

echo "[$(date)] Downloading taxonomy..."
kraken2-build \
  --download-taxonomy \
  --db $KRAKEN2_DB \
  --threads $SLURM_CPUS_PER_TASK

echo "[$(date)] Downloading libraries..."
kraken2-build \
  --download-library bacteria \
  --db $KRAKEN2_DB \
  --threads $SLURM_CPUS_PER_TASK

kraken2-build \
  --download-library archaea \
  --db $KRAKEN2_DB \
  --threads $SLURM_CPUS_PER_TASK

kraken2-build \
  --download-library viral \
  --db $KRAKEN2_DB \
  --threads $SLURM_CPUS_PER_TASK

kraken2-build \
  --download-library human \
  --db $KRAKEN2_DB \
  --threads $SLURM_CPUS_PER_TASK

kraken2-build \
  --download-library UniVec_Core \
  --db $KRAKEN2_DB \
  --threads $SLURM_CPUS_PER_TASK

kraken2-build \
  --download-library plasmid \
  --db $KRAKEN2_DB \
  --threads $SLURM_CPUS_PER_TASK

# Optional: copy in any custom FASTA files
# cp /path/to/custom/*.fna $KRAKEN2_DB/library/

echo "[$(date)] Building the Kraken2 database (this may take a while)..."
kraken2-build \
  --build \
  --db $KRAKEN2_DB \
  --threads $SLURM_CPUS_PER_TASK

echo "[$(date)] Cleaning up intermediate files..."
kraken2-build \
  --clean \
  --db $KRAKEN2_DB

echo "[$(date)] Kraken2 database build complete. Files in $KRAKEN2_DB"