# Sequencing Read Simulator

A powerful tool for simulating Illumina paired-end and Oxford Nanopore (ONT) reads from FASTA files with realistic error profiles and customizable parameters.

## Overview

This tool simulates sequencing reads from a FASTA file containing contigs or genomes. It supports both Illumina paired-end reads and Oxford Nanopore (ONT) long reads. The tool uses BBMap's randomreads.sh to create realistic reads with platform-specific error profiles and customizable parameters.

## Features

- **Multi-Platform Support**: Generate either Illumina paired-end or Oxford Nanopore (ONT) long reads
- **Realistic Read Simulation**: Generates reads with authentic platform-specific error profiles
- **Customizable Parameters**: Adjust insert size, read length, sequencing depth, and error rates
- **Copy Number Control**: Specify different abundances for specific contigs using a TSV file
- **Organized Output**: Structured output directories for reads, statistics, and logs
- **Comprehensive Logging**: Detailed logs of all simulation steps and parameters
- **Summary Reports**: Easy-to-read summaries of simulation results

## Requirements

All dependencies are managed through Pixi:
- BBMap for read simulation and genome fragmentation
- SeqKit for sequence manipulation and statistics
- Samtools for file conversion

## Installation

1. Make sure you have Pixi installed:
   ```bash
   curl -fsSL https://pixi.sh/install.sh | bash
   ```

2. Initialize the Pixi environment:
   ```bash
   cd read_simulator
   pixi install
   ```

## Usage

### Basic Usage

Run the simulation with default parameters:
```bash
pixi run simulate
```

This will:
1. Use `test.fna` as the input file
2. Simulate 150bp paired-end reads at 30x coverage with 500bp insert size
3. Output the results to the `output` directory

### Advanced Usage

Customize the simulation with various parameters:

```bash
pixi run simulate-with-args -- \
  --input my_genome.fna \
  --output my_sample \
  --outdir results \
  --insert_size 400 \
  --insert_stdev 50 \
  --read_length 150 \
  --depth 20 \
  --error_rate 0.01 \
  --threads 8
```

### Available Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Input FASTA file with contigs | `test.fna` |
| `--output` | Output prefix for FASTQ files | `simulated` |
| `--outdir` | Output directory | `output` |
| `--platform` | Sequencing platform (`illumina` or `ont`) | `illumina` |
| `--insert_size` | Mean insert size for Illumina (bp) | `500` |
| `--insert_stdev` | Standard deviation of insert size for Illumina | `50` |
| `--read_length` | Read length (bp) | `150` for Illumina, `2000` for ONT |
| `--depth` | Sequencing depth (x) | `30` |
| `--error_rate` | Base sequencing error rate | `0.01` |
| `--threads` | Number of threads to use | `4` |
| `--interleaved` | Create interleaved reads (for Illumina only) | `false` |

### Clean Up

To remove all output files:
```bash
pixi run clean
```

## Output Structure

The tool creates an organized directory structure using the input file basename:

```
# For standard output (paired-end reads):
sample_simulated/                  # Output directory (<basename>_simulated)
├── reads_sim/                     # Contains FASTQ files
│   ├── sample_R1.fastq.gz         # Forward reads (<basename>_R1.fastq.gz)
│   └── sample_R2.fastq.gz         # Reverse reads (<basename>_R2.fastq.gz)
├── stats/                         # Contains statistics files
│   ├── sample_original_stats.txt  # Input FASTA statistics
│   ├── sample_processed_stats.txt # Processed FASTA statistics (if copy numbers used)
│   └── sample_reads_stats.txt     # Output reads statistics
├── log/                           # Contains log files
│   └── sample_20230504_123456.log # Detailed log of the simulation
└── sample_summary.txt             # Summary of the simulation

# For interleaved output (when using --interleaved):
sample_simulated/                  # Output directory (<basename>_simulated)
├── reads_sim/                     # Contains FASTQ files
│   └── sample.fastq.gz            # Interleaved reads (<basename>.fastq.gz)
├── stats/                         # Contains statistics files
│   ├── sample_original_stats.txt  # Input FASTA statistics
│   ├── sample_processed_stats.txt # Processed FASTA statistics (if copy numbers used)
│   └── sample_reads_stats.txt     # Output reads statistics
├── log/                           # Contains log files
│   └── sample_20230504_123456.log # Detailed log of the simulation
└── sample_summary.txt             # Summary of the simulation
```

## How It Works

The simulation process follows these steps:

1. **Input Analysis**: Analyzes the input contigs to determine their characteristics
2. **Copy Number Processing** (optional):
   - Reads a TSV file with contig IDs and copy numbers
   - Duplicates contigs according to their specified copy numbers
   - Creates a processed FASTA file with the duplicated contigs
3. **Read Simulation**: Uses BBMap's randomreads.sh to:
   - Fragment contigs based on the specified insert size distribution
   - Generate paired-end reads with the specified length
   - Apply realistic error profiles to the reads
4. **Output Organization**: Organizes output files into a structured directory hierarchy
5. **Statistics Generation**: Produces comprehensive statistics about the input and output

### Copy Number Control

You can control the relative abundance of specific contigs in your simulation by providing a TSV file with copy numbers. The file should be named `<basename>/<basename>.tsv` where `<basename>` is the name of your input FASTA file without the extension.

For example, if your input file is `test/test.fna`, the copy number file should be `test/test.tsv`.

The TSV file format is simple:
```
contig_id    copy_number
```

For example:
```
Legionella_pneumophila_fraseri_ATCC_33216    5
Bovine_rhinitis_BB_EU236594                  2
```

This means:
- The contig `Legionella_pneumophila_fraseri_ATCC_33216` will be duplicated to have 5x coverage
- The contig `Bovine_rhinitis_BB_EU236594` will be duplicated to have 2x coverage
- All other contigs will have 1x coverage (default)

The tool will:
1. Read the original FASTA file
2. Create additional copies of the specified contigs
3. Use this processed FASTA file for read simulation

This results in higher read coverage for the duplicated contigs, proportional to their copy numbers.

## Examples

### Simulate Illumina Human Microbiome Reads

```bash
pixi run simulate-with-args -- \
  --input reference_genomes.fna \
  --output human_microbiome \
  --platform illumina \
  --depth 15 \
  --insert_size 350
```

### Simulate Illumina Viral Metagenome with Short Inserts

```bash
pixi run simulate-with-args -- \
  --input viral_contigs.fna \
  --output viral_meta \
  --platform illumina \
  --insert_size 200 \
  --insert_stdev 25 \
  --depth 100
```

### Simulate ONT Long Reads

```bash
pixi run simulate-with-args -- \
  --input bacterial_genome.fna \
  --output bacterial_ont \
  --platform ont \
  --read_length 5000 \
  --depth 20
```

### Simulate Ultra-Long ONT Reads

```bash
pixi run simulate-with-args -- \
  --input assembly.fna \
  --output assembly_ultralong \
  --platform ont \
  --read_length 50000 \
  --depth 5
```

### Simulate Reads with Different Contig Abundances

First, create a TSV file with copy numbers:

```bash
# Create test/test.tsv with copy numbers
echo -e "Legionella_pneumophila_fraseri_ATCC_33216\t5\nBovine_rhinitis_BB_EU236594\t2" > test/test.tsv
```

Then run the simulation:

```bash
pixi run simulate-with-args -- \
  --input test/test.fna \
  --output test_abundance \
  --platform illumina \
  --depth 10
```

This will:
1. Read the copy numbers from `test/test.tsv`
2. Create a processed FASTA file with duplicated contigs
3. Simulate reads with higher coverage for the specified contigs

### Simulate Interleaved Illumina Reads

To generate interleaved reads (where forward and reverse reads are in a single file), use the `--interleaved` flag:

```bash
pixi run simulate-with-args -- \
  --input test/test.fna \
  --output test_interleaved \
  --platform illumina \
  --depth 5 \
  --interleaved
```

This will:
1. Generate paired-end reads from the input FASTA file
2. Create a single interleaved FASTQ file where each forward read is immediately followed by its corresponding reverse read
3. Output the interleaved file to `test_interleaved_simulated/reads_sim/test_interleaved.fastq.gz`

Interleaved reads are useful for tools that expect paired-end reads in a single file.

## Troubleshooting

### Common Issues

1. **Error: Input file does not exist**
   - Make sure the input FASTA file path is correct

2. **Low memory issues**
   - For large genomes, increase available memory or reduce depth

3. **Slow simulation**
   - Use more threads with `--threads` parameter
   - Reduce sequencing depth for faster results
