#!/bin/bash
#SBATCH -p medmem
#SBATCH --job-name=align_slice_downsample
#SBATCH --output=align_slice_downsample.log
#SBATCH --error=align_slice_downsample-%j.err
#SBATCH --time=170:00:00
#SBATCH -c 10
#SBATCH --mem=48G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=clazari@ucsc.edu

set -e

# === LOAD MODULES ===
module load bio/samtools
module load aligners/bwa
module load bio/seqtk/1.3

# === INPUT FILES ===
R1="../fastq_run2/fastq/SRR23782967_1.fastq"
R2="../fastq_run2/fastq/SRR23782967_2.fastq"
REF="../resources/genome.fasta"
GENOME_SIZE=2406167875
READ_LENGTH=150
TARGET_DEPTH=2

# === CONFIG ===
CHUNK_SIZE=40000000 #this is N of my reads paired on the other samples
SEED=100
THREADS=8
BATCH_SIZE=50  # Adjust based on your ulimit

# === WORKING DIRS ===
CHUNK_DIR="fastq_chunks_only"
BAM_DIR="bam_chunks_3"
CHROM_DIR="per_chrom_bams_3"
MERGED_BAM="SRR23782967_downsampled.sorted.bam"
OUTPUT_DIR="$(pwd)/final_output_3"
MERGED_BAM="$OUTPUT_DIR/SRR23782967_downsampled.sorted.bam"

mkdir -p "$CHUNK_DIR" "$BAM_DIR" "$CHROM_DIR" "$OUTPUT_DIR" || {
    echo "[ERROR] Failed to create one or more output directories. Exiting."
    exit 1
}

echo "[DEBUG] Created directories:"
ls -ld "$CHUNK_DIR" "$BAM_DIR" "$CHROM_DIR" "$OUTPUT_DIR"

# === STEP 1: SLICE FASTQ FILES INTO CHUNKS ===
echo "[INFO] Slicing FASTQ files into chunks of $CHUNK_SIZE read pairs..."
seqtk sample -s$SEED "$R1" $CHUNK_SIZE > "$CHUNK_DIR/R1_chunk_01.fastq"
seqtk sample -s$SEED "$R2" $CHUNK_SIZE > "$CHUNK_DIR/R2_chunk_01.fastq"