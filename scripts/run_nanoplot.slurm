#!/bin/bash
#SBATCH --job-name=nanoplot_lizard
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

cd /ibex/user/majnouym/cs249_a2/lizard_asembly

module purge
conda init
conda activate cs249

cd /ibex/user/majnouym/cs249_a2/lizard_asembly

NanoPlot -t 64 -o /ibex/user/majnouym/cs249_a2/lizard_asembly/NanoPlot --fastq   /ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_rev.fastq.gz /ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz