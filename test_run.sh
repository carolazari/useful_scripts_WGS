#!/bin/bash
#SBATCH -p medmem
#SBATCH --job-name=align_slice_downsample
#SBATCH --output=align_slice_downsample.log
#SBATCH --error=align_slice_downsample-%j.err
#SBATCH --time=160:00:00
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
TARGET_DEPTH=8

# === CONFIG ===
CHUNK_SIZE=10000
SEED=100
THREADS=4

# === WORKING DIRS ===
CHUNK_DIR="fastq_chunks_test"
BAM_DIR="bam_chunks_test"
CHROM_DIR="per_chrom_bams_test"
OUTPUT_DIR="final_output_test"
MERGED_BAM="$OUTPUT_DIR/SRR23782967_downsampled.sorted.bam"

mkdir -p "$CHUNK_DIR" "$BAM_DIR" "$CHROM_DIR" "$OUTPUT_DIR"

# === STEP 1: SLICE FASTQ FILES INTO CHUNKS ===
echo "[INFO] Slicing FASTQ files into chunks of $CHUNK_SIZE read pairs..."
seqtk sample -s$SEED "$R1" $CHUNK_SIZE > "$CHUNK_DIR/R1_chunk_01.fastq"
seqtk sample -s$SEED "$R2" $CHUNK_SIZE > "$CHUNK_DIR/R2_chunk_01.fastq"

# === STEP 2: INDEX REFERENCE GENOME ===
echo "[INFO] Indexing reference genome..."
bwa index "$REF"
samtools faidx "$REF"

# === STEP 3: ALIGN, SORT, VERIFY, AND INDEX EACH CHUNK ===
echo "[INFO] Aligning FASTQ chunks..."
for i in "$CHUNK_DIR"/R1_chunk_*.fastq; do
    base=$(basename "$i" | sed 's/R1_chunk_//;s/.fastq//')
    R1_CHUNK="$CHUNK_DIR/R1_chunk_${base}.fastq"
    R2_CHUNK="$CHUNK_DIR/R2_chunk_${base}.fastq"
    OUT_BAM="$BAM_DIR/chunk_${base}.bam"
    SORTED_BAM="$BAM_DIR/chunk_${base}.sorted.bam"

    echo "  Aligning chunk $base..."
    bwa mem -t $THREADS "$REF" "$R1_CHUNK" "$R2_CHUNK" | samtools view -b -o "$OUT_BAM" -

    echo "[INFO] Verifying BAM integrity..."
    if ! samtools quickcheck -v "$OUT_BAM"; then
        echo "[WARNING] BAM invalid, re-running alignment..."
        rm -f "$OUT_BAM"
        bwa mem -t 4 "$REF" "$R1_CHUNK" "$R2_CHUNK" | samtools view -b -o "$OUT_BAM" -
        if ! samtools quickcheck -v "$OUT_BAM"; then
            echo "[ERROR] BAM still invalid. Exiting."
            exit 1
        fi
    fi
    echo "[INFO] $OUT_BAM is valid."

    echo "  Sorting and indexing chunk $base..."
    samtools sort -@ $THREADS -o "$SORTED_BAM" "$OUT_BAM"
    samtools index "$SORTED_BAM"

done

# === STEP 4: CALCULATE COVERAGE AND DOWNSAMPLE ===
echo "[INFO] Calculating coverage and sampling fraction..."
LAST_SORTED_BAM=$(ls -v "$BAM_DIR"/chunk_*.sorted.bam | tail -n 1)

if [ ! -f "$LAST_SORTED_BAM" ]; then
    echo "[ERROR] No BAM files found for coverage calculation."
    exit 1
fi

TOTAL_READS=$(samtools view -c "$LAST_SORTED_BAM")
CURRENT_DEPTH=$(echo "$TOTAL_READS * $READ_LENGTH / $GENOME_SIZE" | bc -l)
FRACTION=$(echo "$TARGET_DEPTH / $CURRENT_DEPTH" | bc -l)
DECIMAL=$(printf "%.3f" "$FRACTION" | cut -d'.' -f2)

echo "  Total reads: $TOTAL_READS"
echo "  Estimated depth: $CURRENT_DEPTH"
echo "  Target depth: $TARGET_DEPTH"
echo "  Downsampling fraction: $FRACTION"

# === STEP 5: PER-CHROMOSOME PROCESSING ===
echo "[INFO] Getting chromosome list..."
readarray -t CHROMOSOMES < <(samtools idxstats "$LAST_SORTED_BAM" | cut -f1 | grep -v '\*')
declare -p CHROMOSOMES

if [ ${#CHROMOSOMES[@]} -eq 0 ]; then
    echo "[ERROR] No chromosomes found. BAM might not be sorted or indexed correctly."
    exit 1
fi

for c in "${CHROMOSOMES[@]}"; do
    echo "Processing $c"
    OUT1="$CHROM_DIR/${c}.aligned.bam"
    OUT2="$CHROM_DIR/${c}.downsampled.bam"
    OUT3="$CHROM_DIR/${c}.downsampled.sorted.bam"

    samtools view -b "$LAST_SORTED_BAM" "$c" > "$OUT1" || echo "[WARN] Could not extract $c"
    samtools view -@ 4 -f 3 -s 42."$DECIMAL" -b "$OUT1" > "$OUT2" || echo "[WARN] Downsampling failed for $c"
    samtools sort -@ 4 -o "$OUT3" "$OUT2" || echo "[WARN] Sorting failed for $c"

    if [ -f "$OUT3" ]; then
        samtools index "$OUT3"
    else
        echo "[WARN] $OUT3 was not created!"
    fi

done

# === STEP 6: MERGE DOWNSAMPLED PER-CHROM BAMs IN BATCHES ===
echo "[INFO] Merging downsampled per-chromosome BAMs in batches..."
TMP_DIR="$CHROM_DIR/tmp_merge_test"
mkdir -p "$TMP_DIR"

FILES=($(ls "$CHROM_DIR"/*.downsampled.sorted.bam 2>/dev/null))
TOTAL_FILES=${#FILES[@]}
BATCH_SIZE=100
MERGED_LIST=()

if [ "$TOTAL_FILES" -eq 0 ]; then
    echo "[ERROR] No BAM files found for merging. Exiting."
    exit 1
fi

HEADER_SOURCE="${FILES[0]}"

for ((i=0; i<TOTAL_FILES; i+=BATCH_SIZE)); do
    BATCH_FILES=("${FILES[@]:i:BATCH_SIZE}")
    BATCH_INDEX=$((i / BATCH_SIZE + 1))
    BATCH_OUT="$TMP_DIR/merge_batch_${BATCH_INDEX}.bam"

    if [ "${#BATCH_FILES[@]}" -eq 0 ]; then
        echo "[WARNING] Empty batch at index $BATCH_INDEX, skipping."
        continue
    fi

    echo "  [INFO] Merging batch $BATCH_INDEX with ${#BATCH_FILES[@]} files..."
    samtools merge -f -@ $THREADS -h "$HEADER_SOURCE" "$BATCH_OUT" "${BATCH_FILES[@]}"
    MERGED_LIST+=("$BATCH_OUT")
done

# Final merge from batches
echo "[INFO] Performing final merge..."
samtools merge -f -@ $THREADS -h "$HEADER_SOURCE" "$MERGED_BAM" "${MERGED_LIST[@]}"
samtools index "$MERGED_BAM"

# === FINAL OUTPUT ===
echo "[INFO] Pipeline complete. Output files:"
echo " - $MERGED_BAM"
echo " - ${MERGED_BAM}.bai"
