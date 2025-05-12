#!/bin/bash
#SBATCH --job-name=qt_1
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH -o logs/%A.out
#SBATCH -e logs/%A.err
#SBATCH --mail-type=ALL

module load fastp

fastp -i /ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz -o /ibex/user/majnouym/cs249_a2/new_assembly/filtered_reads/filtered_seq.fastq.gz --qualified_quality_phred 20 --length_required 1 --disable_length_filtering --thread 16 --html fastp_report_2.html --json  fastp_report_2.json