#!/bin/bash
#SBATCH --job-name=hifiasm_lizard
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

echo "Running on $(hostname)"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"

module purge
module load hifiasm/0.25.0 

cd /ibex/user/majnouym/cs249_a2/hifiasm

hifiasm -o hifiasm_lizard -t 32 /ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_rev.fastq.gz /ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz

echo "End time: $(date)"