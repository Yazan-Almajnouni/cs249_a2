reads=( \
  /ibex/user/majnouym/cs249_a2/new_assembly/filtered_reads/filtered_rev.fastq.gz \
  /ibex/user/majnouym/cs249_a2/new_assembly/filtered_reads/filtered_seq.fastq.gz \
)

outs=( \
  lizard_rev.k21.meryl \
  lizard_seq.k21.meryl \
)

# Loop over indices 0 and 1
for i in "${!reads[@]}"; do
  read="${reads[i]}"
  out="${outs[i]}"
  echo "Submitting job for read=$read → out=$out"
  sbatch /ibex/user/majnouym/cs249_a2/scripts/run_myrel_2.slurm /ibex/user/majnouym/cs249_a2/new_assembly/meryl "$out" "$read" 
done