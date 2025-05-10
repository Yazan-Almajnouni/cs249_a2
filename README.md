
Tasks 1.1-1.3

LLM: GPT-o4-mini
Prompt: Lets do the following tasks in steps. In the first step, only do Task 1.1. Keep the other tasks in mind though for things like output format and such things.
Tasks: [tasks inserted here one by one]

LLM output is used with slight modifications.


# De Bruijn Graph Assembler

A simple De Bruijn Graph‐based genome assembler.

## Requirements

- Python 3.6 or higher

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

Given `sample.fastq`:

Running

python src/dbg_assembler.py \
  -i b.fastq \
  -k 40 \
  --min-count 2 \
  --gfa b.k40.graph.gfa \
  -o b.k40.contigs.fasta

might produce `b.k40.contigs.fasta` and `b.k40.graph.gfa`


sources:
- https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf
