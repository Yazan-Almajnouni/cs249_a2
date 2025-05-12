#!/usr/bin/env bash

set -euo pipefail

module load quast

# List of assemblies
assemblies=(
  /ibex/user/majnouym/cs249_a2/outputs/task_1_3_3/hiseq/error/assymbly_dbg.fasta
  /ibex/user/majnouym/cs249_a2/outputs/task_1_3_3/hiseq/no_error/assymbly_dbg.fasta
  /ibex/user/majnouym/cs249_a2/outputs/task_1_3_3/ont/error/assymbly_dbg.fasta
  /ibex/user/majnouym/cs249_a2/outputs/task_1_3_3/ont/no_error/assymbly_dbg.fasta
  /ibex/user/majnouym/cs249_a2/outputs/task_1_3_3/hiseq/error/assymbly_olc.fasta
  /ibex/user/majnouym/cs249_a2/outputs/task_1_3_3/hiseq/no_error/assymbly_olc.fasta
  /ibex/user/majnouym/cs249_a2/outputs/task_1_3_3/ont/error/assymbly_olc.fasta
  /ibex/user/majnouym/cs249_a2/outputs/task_1_3_3/ont/no_error/assymbly_olc.fasta
)

# Base output directory for all QUAST runs
base_outdir="/ibex/user/majnouym/cs249_a2/outputs/task_1_3_3/quast_reports"

for asm in "${assemblies[@]}"; do
  # strip everything up to and including "task_1_3_3/"
  rel_path="${asm#*/task_1_3_3/}"
  # get the directory part (e.g. "hiseq/error")
  dir_part="$(dirname "$rel_path")"
  # get the filename without extension (e.g. "assymbly_dbg")
  base_name="$(basename "$rel_path" .fasta)"

  # construct the output directory
  outdir="$base_outdir/$dir_part/$base_name"

  echo "Running QUAST on $asm -> $outdir"
  mkdir -p "$outdir"
  quast.py "$asm" -o "$outdir" -r /ibex/user/majnouym/cs249_a2/data/synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna --min-contig 0
done