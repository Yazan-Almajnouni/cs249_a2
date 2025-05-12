# Source Code For CS249 Assignment 2

## Requirements

- Python 3.9 or higher

## De Bruijn Graph Assembler

A simple De Bruijn Graph‐based genome assembler.

### Usage

`python src/dbg_assembler.py args`

Arguments:
  - `-i, --input`  One or more FASTQ files
  - `-k, --kmer`   k-mer size 
  - `--min-count`  Minimum k-mer count to include in graph (default: 1)
  - `-o, --output` Output FASTA file for contigs
  - `--gfa <file>`     If set, also output the raw de Bruijn graph in GFA1
                   (for Bandage, etc.)

### Example

Running

```
python src/dbg_assembler.py \
  -i data/toy_dataset/reads_b.fastq \
  -k 40 \
  --gfa outputs/task_1_3_1/reads_b.gfa \
  -o outputs/task_1_3_1/assymbly_b.fasta
```

produces `outputs/task_1_3_1/assymbly_b.fasta` and `outputs/task_1_3_1/reads_b.gfa`

## Overlap-Layout-Consensus Assembler

A simple Overlap-Layout-Consensus-based genome assembler.

### Usage

`python src/olc_assembler.py args`

Arguments:
  - `-i, --input`        One or more FASTQ files
  - `-n, --min-overlap`  Minimum overlap length in the overlap graph
  - `-o, --output`       Output FASTA file for contigs

### Example

Running

```
python src/olc_assembler.py \
  -i data/toy_dataset/reads_b.fastq \
  -n 50
  -o outputs/misc/assymbly_b_olc.fasta
```

produces `outputs/misc/assymbly_b_olc.fasta`

## data

Contains the data for task 1.

## scripts

Contains scripts to reproduce the analysis and evaluations for the assignment. The scripts run with no command line arguments unless stated otherwise. 
Some of the scripts are specific to my ibex folder since they use absolute path.

- `build_k8.slurm`: used to build kraken2 standard-8 database 
- `build_krakendb.slurm`: used to build kraken2 standard database
- `hifiasm_to_fasta`: used to convert hifiasm output contigs to FASTA format.
- `merge_meryl.slurm`: used to merge meryl databases from the two PacBio reads files.
- `merge_meryl.slurm`: same as above, but for the filtered reads.
- `plot_k-mer.py`: plots the k-mer distribution of the first assembly for task 2.2.
- `run_busco.slurm`: runs busco on the assembly to report gene completeness.
- `run_fastp_1.slurm`: filters the `lizard_liver_rev.fastq.gz` file for reads below 20 mean quality.
- `run_fastp_2.slurm`: same as above, but for `lizard_liver_seq.fastq.gz` file.
- `run_hifiasm_2.slurm`: assembles the lizard from the filtered reads.
- `run_hifiasm.slurm`: assembles the lizard from the old reads.
- `run_inspector.slurm`: runs inspector on the assemby.
- `run_kraken.slurm`: runs kraken2 to classify reads in the original files.
- `run_merqury.slurm`: calculates QV score for first assembly.
- `run_merqury2.slurm` calculates QV score for new assembly.
- `run_meryl_new_assembly.sh`: calls `run_meryl_2.slurm` on the filtered reads files.
- `run_meryl.sh`: calls `run_meryl_2.slurm` on the original reads files.
- `run_meryl_2.slurm`: builds meryl database. takes 3 args: dir, out, read. dir is the directory where the reads are, out is the name of the output database, and read is the reads file name.
- `run_nanoplot.slurm`: runs nanoplot on the original reads to perform quality assessment.
- `run_quast_1_3_3.sh`: runs quast on the dbg and olc assemblies for task 1.3.3.
- `run_quast_1_3_4.sh`: runs quast on the SPAdes assemblies for task 1.3.4.
- `run_spades_error`: runs spades on the ONT reads with error.
- `run_spades_no_error`: runs spades on the ONT reads with no error.



