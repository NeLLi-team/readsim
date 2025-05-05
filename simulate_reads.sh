#!/bin/bash
set -e

# Default parameters
INPUT_FILE="test.fna"
OUTPUT_PREFIX=""  # Will be set based on input file basename if not provided
OUTPUT_DIR=""     # Will be set to OUTPUT_PREFIX if not provided
PLATFORM="illumina"  # Options: illumina, ont
INSERT_SIZE=500
INSERT_STDEV=50
READ_LENGTH=""  # Will be set based on platform
DEPTH=30
ERROR_RATE=0.01
THREADS=4
INTERLEAVED=false  # Whether to create interleaved reads (for Illumina only)

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --input)
            INPUT_FILE="$2"
            shift 2
            ;;
        --output)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        --outdir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --platform)
            PLATFORM=$(echo "$2" | tr '[:upper:]' '[:lower:]')  # Convert to lowercase
            if [[ "$PLATFORM" != "illumina" && "$PLATFORM" != "ont" ]]; then
                echo "Error: Platform must be either 'illumina' or 'ont'"
                exit 1
            fi
            shift 2
            ;;
        --insert_size)
            INSERT_SIZE="$2"
            shift 2
            ;;
        --insert_stdev)
            INSERT_STDEV="$2"
            shift 2
            ;;
        --read_length)
            READ_LENGTH="$2"
            shift 2
            ;;
        --depth)
            DEPTH="$2"
            shift 2
            ;;
        --error_rate)
            ERROR_RATE="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --interleaved)
            INTERLEAVED=true
            shift 1
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --input FILE       Input FASTA file with contigs (default: test.fna)"
            echo "  --output PREFIX    Output prefix for FASTQ files (default: simulated)"
            echo "  --outdir DIR       Output directory (default: output)"
            echo "  --platform STR     Sequencing platform: 'illumina' or 'ont' (default: illumina)"
            echo "  --insert_size INT  Mean insert size for Illumina (default: 500)"
            echo "  --insert_stdev INT Standard deviation of insert size for Illumina (default: 50)"
            echo "  --read_length INT  Read length (default: 150 for Illumina, 2000 for ONT)"
            echo "  --depth FLOAT      Sequencing depth (default: 30)"
            echo "  --error_rate FLOAT Sequencing error rate (default: 0.01)"
            echo "  --threads INT      Number of threads to use (default: 4)"
            echo "  --interleaved      Create interleaved reads (for Illumina only)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE does not exist"
    exit 1
fi

# Extract basename without extension
BASENAME=$(basename "$INPUT_FILE" | sed 's/\.[^.]*$//')

# Check if a copy number TSV file exists
INPUT_DIR=$(dirname "$INPUT_FILE")
TSV_FILE="$INPUT_DIR/$BASENAME.tsv"

# Create temporary directory early
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

PROCESSED_FASTA="$TMP_DIR/${BASENAME}_processed.fna"

# Function to process FASTA file with copy numbers
process_fasta_with_copy_numbers() {
    local input_fasta=$1
    local tsv_file=$2
    local output_fasta=$3

    # If TSV file doesn't exist, just copy the input file
    if [ ! -f "$tsv_file" ]; then
        echo "No copy number file found at $tsv_file. Using original FASTA file."
        cp "$input_fasta" "$output_fasta"
        return
    fi

    echo "Found copy number file: $tsv_file"
    echo "Processing FASTA file with copy numbers..."

    # Copy the original FASTA file first
    cp "$input_fasta" "$output_fasta"

    # Read the TSV file and duplicate contigs as needed
    while IFS=$'\t' read -r contig_id copy_number || [ -n "$contig_id" ]; do
        # Skip empty lines or comments
        if [[ -z "$contig_id" || "$contig_id" == \#* ]]; then
            continue
        fi

        # Calculate how many additional copies we need (original already exists once)
        additional_copies=$((copy_number - 1))

        if [ "$additional_copies" -gt 0 ]; then
            echo "  Adding $additional_copies additional copies of contig: $contig_id"

            # Extract the contig sequence
            contig_seq=$(awk -v RS=">" -v contig="$contig_id" '$1 == contig {print ">"$0}' "$input_fasta")

            if [ -z "$contig_seq" ]; then
                echo "  Warning: Contig $contig_id not found in FASTA file"
                continue
            fi

            # Add the additional copies
            for i in $(seq 1 $additional_copies); do
                # Modify the header to indicate it's a copy
                modified_header=$(echo "$contig_seq" | head -n 1 | sed "s/>/>${contig_id}_copy$i /")
                modified_seq=$(echo "$contig_seq" | tail -n +2)
                echo "$modified_header" >> "$output_fasta"
                echo "$modified_seq" >> "$output_fasta"
            done
        fi
    done < "$tsv_file"

    # Count the number of sequences in the original and processed files
    orig_count=$(grep -c "^>" "$input_fasta")
    proc_count=$(grep -c "^>" "$output_fasta")
    added_count=$((proc_count - orig_count))

    echo "FASTA processing complete: Added $added_count copies of contigs"
    echo "Original sequences: $orig_count, Processed sequences: $proc_count"
}

# Set OUTPUT_PREFIX based on input file basename if not provided
if [ -z "$OUTPUT_PREFIX" ]; then
    OUTPUT_PREFIX="$BASENAME"
    echo "Output prefix not specified, using input file basename: $OUTPUT_PREFIX"
fi

# Set OUTPUT_DIR to <basename>_simulated if not provided
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="${BASENAME}_simulated"
    echo "Output directory not specified, using: $OUTPUT_DIR"
fi

# If output prefix is specified but output dir is not, create a custom output dir
if [ "$OUTPUT_PREFIX" != "$BASENAME" ] && [ "$OUTPUT_DIR" = "${BASENAME}_simulated" ]; then
    OUTPUT_DIR="${OUTPUT_PREFIX}_simulated"
    echo "Custom output prefix detected, using directory: $OUTPUT_DIR"
fi

# Set read length based on platform if not provided
if [ -z "$READ_LENGTH" ]; then
    if [ "$PLATFORM" = "illumina" ]; then
        READ_LENGTH=150
        echo "Read length not specified, using default for Illumina: $READ_LENGTH bp"
    else  # ONT
        READ_LENGTH=2000
        echo "Read length not specified, using default for ONT: $READ_LENGTH bp"
    fi
fi

# Create output directory structure
mkdir -p "$OUTPUT_DIR/reads_sim"
mkdir -p "$OUTPUT_DIR/stats"
mkdir -p "$OUTPUT_DIR/log"

# Get absolute paths
INPUT_FILE_ABS=$(realpath "$INPUT_FILE")
OUTPUT_DIR_ABS=$(realpath "$OUTPUT_DIR")

# Debug paths (commented out)
# echo "Debug paths:"
# echo "TMP_DIR: $TMP_DIR"
# echo "INPUT_FILE: $INPUT_FILE"
# echo "TSV_FILE: $TSV_FILE"
# echo "PROCESSED_FASTA: $PROCESSED_FASTA"

# Process the FASTA file with copy numbers if TSV file exists
process_fasta_with_copy_numbers "$INPUT_FILE" "$TSV_FILE" "$PROCESSED_FASTA"

# Use the processed FASTA file for read simulation
SIMULATION_INPUT="$PROCESSED_FASTA"

# Generate timestamp for log files
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$OUTPUT_DIR/log/${OUTPUT_PREFIX}_${TIMESTAMP}.log"

# Log simulation parameters
{
    echo "=== Read Simulation Parameters ==="
    echo "Timestamp: $(date)"
    echo "Input file: $INPUT_FILE_ABS"
    echo "Output directory: $OUTPUT_DIR_ABS"
    echo "Output prefix: $OUTPUT_PREFIX"
    echo "Platform: $PLATFORM"
    if [ "$PLATFORM" = "illumina" ]; then
        echo "Insert size: $INSERT_SIZE ± $INSERT_STDEV"
        echo "Interleaved reads: $INTERLEAVED"
    fi
    echo "Read length: $READ_LENGTH"
    echo "Sequencing depth: ${DEPTH}x"
    echo "Error rate: $ERROR_RATE"
    echo "Threads: $THREADS"
    echo "Temporary directory: $TMP_DIR"
    echo "===================================="
} | tee "$LOG_FILE"

# Step 1: Analyze input contigs
echo -e "\n# Step 1: Analyzing input contigs..." | tee -a "$LOG_FILE"
seqkit stats "$INPUT_FILE" > "$OUTPUT_DIR/stats/${OUTPUT_PREFIX}_original_stats.txt"
cat "$OUTPUT_DIR/stats/${OUTPUT_PREFIX}_original_stats.txt" | tee -a "$LOG_FILE"

# If we processed the FASTA file with copy numbers, analyze that too
if [ -f "$TSV_FILE" ]; then
    echo -e "\n# Step 1b: Analyzing processed contigs with copy numbers..." | tee -a "$LOG_FILE"
    seqkit stats "$SIMULATION_INPUT" > "$OUTPUT_DIR/stats/${OUTPUT_PREFIX}_processed_stats.txt"
    cat "$OUTPUT_DIR/stats/${OUTPUT_PREFIX}_processed_stats.txt" | tee -a "$LOG_FILE"
fi

# Step 2: Generate reads based on platform
echo -e "\n# Step 2: Generating $PLATFORM reads..." | tee -a "$LOG_FILE"

if [ "$PLATFORM" = "illumina" ]; then
    # Use randomreads.sh to generate simulated Illumina paired-end reads
    {
        randomreads.sh \
            ref="$SIMULATION_INPUT" \
            out="$TMP_DIR/${OUTPUT_PREFIX}1.fq" \
            out2="$TMP_DIR/${OUTPUT_PREFIX}2.fq" \
            paired=t \
            length="$READ_LENGTH" \
            mininsert="$((INSERT_SIZE - 2*INSERT_STDEV))" \
            maxinsert="$((INSERT_SIZE + 2*INSERT_STDEV))" \
            coverage="$DEPTH" \
            adderrors=t \
            maxsnps=0 \
            maxinss=0 \
            maxdels=0 \
            maxsubs=0 \
            snprate=0.001 \
            insrate=0.0001 \
            delrate=0.0001 \
            subrate=0.0001 \
            seed="$RANDOM" \
            illuminanames=t \
            addpairnum=t
    } 2>&1 | tee -a "$LOG_FILE"
else
    # Use randomreads.sh to generate simulated ONT long reads
    {
        randomreads.sh \
            ref="$SIMULATION_INPUT" \
            out="$TMP_DIR/${OUTPUT_PREFIX}1.fq" \
            paired=f \
            length="$READ_LENGTH" \
            coverage="$DEPTH" \
            adderrors=t \
            maxsnps=0 \
            maxinss=0 \
            maxdels=0 \
            maxsubs=0 \
            snprate=0.05 \
            insrate=0.01 \
            delrate=0.01 \
            subrate=0.01 \
            seed="$RANDOM"
    } 2>&1 | tee -a "$LOG_FILE"
fi

echo -e "\nSimulated reads generated successfully." | tee -a "$LOG_FILE"

# Step 3: Compress output files
echo -e "\n# Step 3: Compressing output files..." | tee -a "$LOG_FILE"

if [ "$PLATFORM" = "illumina" ] && [ "$INTERLEAVED" = true ]; then
    # Create interleaved FASTQ file for Illumina paired-end reads
    echo -e "Creating interleaved FASTQ file..." | tee -a "$LOG_FILE"

    # Create a temporary interleaved file by alternating R1 and R2 reads
    paste <(cat "$TMP_DIR/${OUTPUT_PREFIX}1.fq" | paste - - - -) <(cat "$TMP_DIR/${OUTPUT_PREFIX}2.fq" | paste - - - -) | \
    tr '\t' '\n' > "$TMP_DIR/${OUTPUT_PREFIX}.fq"

    # Compress the interleaved file
    gzip -c "$TMP_DIR/${OUTPUT_PREFIX}.fq" > "$OUTPUT_DIR/reads_sim/${OUTPUT_PREFIX}.fastq.gz"

    # Step 4: Generate statistics for simulated reads
    echo -e "\n# Step 4: Generating statistics for simulated reads..." | tee -a "$LOG_FILE"
    seqkit stats "$OUTPUT_DIR/reads_sim/${OUTPUT_PREFIX}.fastq.gz" > "$OUTPUT_DIR/stats/${OUTPUT_PREFIX}_reads_stats.txt"
elif [ "$PLATFORM" = "illumina" ]; then
    # Standard Illumina paired-end reads (non-interleaved)
    gzip -c "$TMP_DIR/${OUTPUT_PREFIX}1.fq" > "$OUTPUT_DIR/reads_sim/${OUTPUT_PREFIX}_R1.fastq.gz"
    gzip -c "$TMP_DIR/${OUTPUT_PREFIX}2.fq" > "$OUTPUT_DIR/reads_sim/${OUTPUT_PREFIX}_R2.fastq.gz"

    # Step 4: Generate statistics for simulated reads
    echo -e "\n# Step 4: Generating statistics for simulated reads..." | tee -a "$LOG_FILE"
    seqkit stats "$OUTPUT_DIR/reads_sim/${OUTPUT_PREFIX}_R1.fastq.gz" "$OUTPUT_DIR/reads_sim/${OUTPUT_PREFIX}_R2.fastq.gz" > "$OUTPUT_DIR/stats/${OUTPUT_PREFIX}_reads_stats.txt"
else
    # ONT single-end reads
    gzip -c "$TMP_DIR/${OUTPUT_PREFIX}1.fq" > "$OUTPUT_DIR/reads_sim/${OUTPUT_PREFIX}_R1.fastq.gz"

    # Step 4: Generate statistics for simulated reads (ONT single-end)
    echo -e "\n# Step 4: Generating statistics for simulated reads..." | tee -a "$LOG_FILE"
    seqkit stats "$OUTPUT_DIR/reads_sim/${OUTPUT_PREFIX}_R1.fastq.gz" > "$OUTPUT_DIR/stats/${OUTPUT_PREFIX}_reads_stats.txt"
fi

cat "$OUTPUT_DIR/stats/${OUTPUT_PREFIX}_reads_stats.txt" | tee -a "$LOG_FILE"

# Create a summary file with key information
SUMMARY_FILE="$OUTPUT_DIR/${OUTPUT_PREFIX}_summary.txt"
{
    echo "# Read Simulation Summary"
    echo "Timestamp: $(date)"
    echo "Input file: $INPUT_FILE_ABS"

    # Add information about copy number processing if applicable
    if [ -f "$TSV_FILE" ]; then
        echo "Copy number file: $TSV_FILE"
        echo "Processed FASTA: $PROCESSED_FASTA"

        # Count original and processed sequences
        orig_count=$(grep -c "^>" "$INPUT_FILE")
        proc_count=$(grep -c "^>" "$PROCESSED_FASTA")
        added_count=$((proc_count - orig_count))

        echo "Copy number processing: Added $added_count copies of contigs"
        echo "  - Original sequences: $orig_count"
        echo "  - Processed sequences: $proc_count"
    fi

    echo "Platform: $PLATFORM"
    echo "Parameters:"
    if [ "$PLATFORM" = "illumina" ]; then
        echo "  - Insert size: $INSERT_SIZE ± $INSERT_STDEV bp"
    fi
    echo "  - Read length: $READ_LENGTH bp"
    echo "  - Sequencing depth: ${DEPTH}x"
    echo "  - Error rate: $ERROR_RATE"
    echo ""
    echo "Output files:"
    if [ "$PLATFORM" = "illumina" ] && [ "$INTERLEAVED" = true ]; then
        echo "  - Interleaved reads: $OUTPUT_DIR_ABS/reads_sim/${OUTPUT_PREFIX}.fastq.gz"
    else
        echo "  - Forward reads: $OUTPUT_DIR_ABS/reads_sim/${OUTPUT_PREFIX}_R1.fastq.gz"
        if [ "$PLATFORM" = "illumina" ]; then
            echo "  - Reverse reads: $OUTPUT_DIR_ABS/reads_sim/${OUTPUT_PREFIX}_R2.fastq.gz"
        fi
    fi
    echo "  - Statistics: $OUTPUT_DIR_ABS/stats/${OUTPUT_PREFIX}_reads_stats.txt"
    if [ -f "$TSV_FILE" ]; then
        echo "  - Original FASTA stats: $OUTPUT_DIR_ABS/stats/${OUTPUT_PREFIX}_original_stats.txt"
        echo "  - Processed FASTA stats: $OUTPUT_DIR_ABS/stats/${OUTPUT_PREFIX}_processed_stats.txt"
    fi
    echo "  - Log file: $OUTPUT_DIR_ABS/log/${OUTPUT_PREFIX}_${TIMESTAMP}.log"
    echo ""
    echo "Read statistics:"
    cat "$OUTPUT_DIR/stats/${OUTPUT_PREFIX}_reads_stats.txt"
} > "$SUMMARY_FILE"

echo -e "\nSimulation complete!" | tee -a "$LOG_FILE"
echo -e "\nOutput structure:" | tee -a "$LOG_FILE"
if [ "$PLATFORM" = "illumina" ] && [ "$INTERLEAVED" = true ]; then
    echo "  $OUTPUT_DIR/reads_sim/${OUTPUT_PREFIX}.fastq.gz (Interleaved reads)" | tee -a "$LOG_FILE"
else
    echo "  $OUTPUT_DIR/reads_sim/${OUTPUT_PREFIX}_R1.fastq.gz (Forward reads)" | tee -a "$LOG_FILE"
    if [ "$PLATFORM" = "illumina" ]; then
        echo "  $OUTPUT_DIR/reads_sim/${OUTPUT_PREFIX}_R2.fastq.gz (Reverse reads)" | tee -a "$LOG_FILE"
    fi
fi
if [ -f "$TSV_FILE" ]; then
    echo "  $OUTPUT_DIR/stats/${OUTPUT_PREFIX}_original_stats.txt (Original FASTA statistics)" | tee -a "$LOG_FILE"
    echo "  $OUTPUT_DIR/stats/${OUTPUT_PREFIX}_processed_stats.txt (Processed FASTA statistics)" | tee -a "$LOG_FILE"
else
    echo "  $OUTPUT_DIR/stats/${OUTPUT_PREFIX}_original_stats.txt (Input statistics)" | tee -a "$LOG_FILE"
fi
echo "  $OUTPUT_DIR/stats/${OUTPUT_PREFIX}_reads_stats.txt (Read statistics)" | tee -a "$LOG_FILE"
echo "  $OUTPUT_DIR/log/${OUTPUT_PREFIX}_${TIMESTAMP}.log (Log file)" | tee -a "$LOG_FILE"
echo "  $OUTPUT_DIR/${OUTPUT_PREFIX}_summary.txt (Summary file)" | tee -a "$LOG_FILE"

echo -e "\nTo use these reads in your pipeline, reference them at:"
if [ "$PLATFORM" = "illumina" ] && [ "$INTERLEAVED" = true ]; then
    echo "  $OUTPUT_DIR_ABS/reads_sim/${OUTPUT_PREFIX}.fastq.gz (interleaved format)"
else
    echo "  $OUTPUT_DIR_ABS/reads_sim/${OUTPUT_PREFIX}_R1.fastq.gz"
    if [ "$PLATFORM" = "illumina" ]; then
        echo "  $OUTPUT_DIR_ABS/reads_sim/${OUTPUT_PREFIX}_R2.fastq.gz"
    else
        echo "  (ONT reads are single-end only)"
    fi
fi
