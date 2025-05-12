#!/bin/bash
#SBATCH --job-name=hifiasm_lizard
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=400G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

echo "Running on $(hostname)"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"

module purge
module load hifiasm/0.25.0 

cd /ibex/user/majnouym/cs249_a2/new_assembly/filtered_reads/hifiasm

hifiasm -o hifiasm_lizard -t 128 /ibex/user/majnouym/cs249_a2/new_assembly/filtered_reads/filtered_rev.fastq.gz /ibex/user/majnouym/cs249_a2/new_assembly/filtered_reads/filtered_seq.fastq.gz

echo "End time: $(date)"