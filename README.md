
Task 1.1

LLM: GPT-o4-mini
Prompt: Lets do the following tasks in steps. In the first step, only do Task 1.1. Keep the other tasks in mind though for things like output format and such things.
Tasks: [task 1.1 inserted here]

LLM output is used almost verbatum with slight modifications.
The LLM was also used to generate parts of the README, specifically the usage and example sections, and the output was slightly modified.


# De Bruijn Graph Assembler

A simple De Bruijn Graph‐based genome assembler.

## Requirements

- Python 3.9 or higher

## Usage

python src/dbg_assembler.py args

Arguments:
  - `-i, --input`  One or more FASTQ files
  - `-k, --kmer`   k-mer size (integer, e.g. 31)
  - `--min-count`  Minimum k-mer count to include in graph (default: 1)
  - `-o, --output` Output FASTA file for contigs
  --gfa <file>     If set, also output the raw de Bruijn graph in GFA1
                   (for Bandage, etc.)

## Example

Running

python src/dbg_assembler.py \
  -i data/toy_dataset/reads_b.fastq \
  -k 35 \
  --gfa outputs/task_1_3_1/reads_b.gfa \
  -o outputs/task_1_3_1/assymbly_b.fasta

produces outputs/task_1_3_1/assymbly_b.fasta and outputs/task_1_3_1/reads_b.gfa


sources:
- https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf

Task 1.2

LLM: GPT-o4-mini

prompt: Implement an OLC assembler that uses suffix trees as shown in these pictures.
The pictures were screenshots from the suffix tree diagrams in https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_olc.pdf
The output was tested and refined for correctness/performance using prompts similar to "The code has issue [issue], identify why this happens and fix it"

# Overlap-Layout-Consensus Assembler

A simple Overlap-Layout-Consensus-based genome assembler.

## Usage

python src/dbg_assembler.py args

Arguments:
  - `-i, --input`        One or more FASTQ files
  - `-n, --min-overlap`  Minimum overlap length in the overlap graph
  - `-o, --output`       Output FASTA file for contigs

## Example

Running

python src/dbg_assembler.py \
  -i data/toy_dataset/reads_b.fastq \
  -n 50
  -o outputs/task_1_3_1/assymbly_b.fasta

produces 

sources:
- https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_olc.pdf

