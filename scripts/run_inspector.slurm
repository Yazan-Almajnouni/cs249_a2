#!/bin/bash
#SBATCH --job-name=inspector_lizard
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=400G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

module load inspector 

inspector.py -c /ibex/user/majnouym/cs249_a2/lizard_asembly/hifiasm/lizard_asm.fasta -r /ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_rev.fastq.gz   /ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz -o /ibex/user/majnouym/cs249_a2/lizard_asembly/inspector_2 -t 128 