#!/usr/bin/env bash

set -euo pipefail

module load quast

cd /ibex/user/majnouym/cs249_a2

# List of assemblies
assemblies=(
  outputs/task_1_3_4/hiseq/error/contigs.fasta
  outputs/task_1_3_4/hiseq/no_error/contigs.fasta
  outputs/task_1_3_4/ont/error/contigs.fasta
  outputs/task_1_3_4/ont/no_error/contigs.fasta
)

# Base output directory for all QUAST runs
base_outdir="quast_reports"

mkdir -p "$base_outdir"

for asm in "${assemblies[@]}"; do
  if [[ ! -f "$asm" ]]; then
    echo "Warning: $asm not found, skipping."
    continue
  fi

  # Extract platform and error-status from path:
  # e.g. outputs/task_1_3_4/hiseq/error/contigs.fasta
  #      -> platform=hiseq, status=error
  platform=$(echo "$asm" | cut -d'/' -f3)
  status=$(echo "$asm"   | cut -d'/' -f4)

  # Prepare output dir for this run
  outdir="${base_outdir}/quast_${platform}_${status}"
  mkdir -p "$outdir"

  echo "Running QUAST on $asm -> $outdir"
  quast.py "$asm" -o "$outdir" -r data/synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna
done

echo "All QUAST runs completed. Reports are in $base_outdir/"
