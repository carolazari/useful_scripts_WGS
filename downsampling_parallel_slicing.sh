#!/bin/bash
#SBATCH -p medmem
#SBATCH --job-name=align_slice_downsample
#SBATCH --output=align_slice_downsample.log
#SBATCH --error=align_slice_downsample-%j.err
#SBATCH --time=160:00:00
#SBATCH -c 10
#SBATCH --mem=16G
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
GENOME_SIZE=2406167875      # e.g. 1 Gb (adjust to your genome)
READ_LENGTH=150             # e.g. 150 bp
TARGET_DEPTH=8              #change to whatever you need

# === CONFIG ===
CHUNK_SIZE=10000        # 1 million read pairs per chunk. default--> Subsample 10000 read pairs from two large paired FASTQ files (remember to use the same random seed to keep pairing):
SEED=100                  #100 this is the default number from the https://github.com/lh3/seqtk . 
THREADS=4

# === WORKING DIRS ===
CHUNK_DIR="fastq_chunks"
BAM_DIR="bam_chunks"
MERGED_BAM="aligned_merged.bam"
mkdir -p "$CHUNK_DIR" "$BAM_DIR". #-p avoids error if the directory already exists

# === STEP 1: SLICE FASTQ FILES INTO CHUNKS ===
echo "[INFO] Slicing FASTQ files into chunks of $CHUNK_SIZE read pairs..."

# seqtk outputs uncompressed FASTQ
seqtk sample -s$SEED "$R1" $CHUNK_SIZE > "$CHUNK_DIR/R1_chunk_01.fastq"
seqtk sample -s$SEED "$R2" $CHUNK_SIZE > "$CHUNK_DIR/R2_chunk_01.fastq"

# You can loop here for multiple chunks if needed (e.g., random or sequential slicing)

# === STEP 2: INDEX REFERENCE GENOME ===
echo "[INFO] Indexing reference genome..."
bwa index "$REF"
samtools faidx "$REF"

# === STEP 3: ALIGN EACH CHUNK AND OUTPUT BAM ===
echo "[INFO] Aligning FASTQ chunks..."
for i in "$CHUNK_DIR"/R1_chunk_*.fastq; do
    base=$(basename "$i" | sed 's/R1_chunk_//;s/.fastq//')
    R1_CHUNK="$CHUNK_DIR/R1_chunk_${base}.fastq"
    R2_CHUNK="$CHUNK_DIR/R2_chunk_${base}.fastq"
    OUT_BAM="$BAM_DIR/chunk_${base}.bam"

    echo "  Aligning chunk $base..."
    bwa mem -t $THREADS "$REF" "$R1_CHUNK" "$R2_CHUNK" | samtools view -b -o "$OUT_BAM" -
done

#=== STEP 3b:VERIFY BAM INTEGRITY ===
echo "[INFO] Verifying BAM integrity..."
if ! samtools quickcheck -v "$OUT_BAM"; then
    echo "[WARNING] BAM invalid, re-running alignment..."
    rm -f "$OUT_BAM"
    bwa mem -t 4 $REF $R1 $R2 | samtools view -b -o "$OUT_BAM" -
    if ! samtools quickcheck -v "$OUT_BAM"; then
        echo "[ERROR] BAM still invalid. Exiting."
        exit 1
    fi
fi
echo "[INFO] "$OUT_BAM" is valid."

# === STEP 4: MERGE ALL CHUNK BAM FILES ===
#echo "[INFO] Merging all chunk BAMs..."
#samtools merge -@ $THREADS "$MERGED_BAM" "$BAM_DIR"/chunk_*.bam
#samtools index "$MERGED_BAM"

# === STEP 4: CALCULATE COVERAGE AND DOWNSAMPLE ===
echo "[INFO] Calculating coverage and sampling fraction..."

TOTAL_READS=$(samtools view -c "$MERGED_BAM")
CURRENT_DEPTH=$(echo "$TOTAL_READS * $READ_LENGTH / $GENOME_SIZE" | bc -l)
FRACTION=$(echo "$TARGET_DEPTH / $CURRENT_DEPTH" | bc -l)

# Extract only the decimal part for samtools -s argument
DECIMAL=$(printf "%.3f" "$FRACTION" | cut -d'.' -f2)

echo "  Total reads: $TOTAL_READS"
echo "  Estimated depth: $CURRENT_DEPTH"
echo "  Target depth: $TARGET_DEPTH"
echo "  Downsampling fraction: $FRACTION"


# === STEP 5: DOWNSAMPLE BAMS ===
#echo "[INFO] Downsampling BAM..."
#samtools view -@ $THREADS -f 3 -s 42."$DECIMAL" -b "$MERGED_BAM" > downsampled.bam
#samtools sort -@ $THREADS -o downsampled.sorted.bam downsampled.bam
#samtools index downsampled.sorted.bam

# === STEP 6: PER-CHROMOSOME PROCESSING ===

# Get chromosome names into a Bash array
echo "[INFO] Getting chromosome list..."
readarray -t CHROMOSOMES < <(samtools idxstats "$OUT_BAM" | cut -f1 | grep -v '\*')

declare -p CHROMOSOMES

CHROM_DIR="per_chrom_bams"
mkdir -p "$CHROM_DIR" #-p avoids error if the directory already exists

#for c in "${CHROMOSOMES[@]}"; do
#    echo "[INFO] Processing chromosome $c"
#    samtools view -b "$OUT_BAM" "$c" > "$CHROM_DIR/${c}.bam"
#    samtools sort -@ $THREADS -o "$CHROM_DIR/${c}.sorted.bam" "$CHROM_DIR/${c}.bam"
#    samtools index "$CHROM_DIR/${c}.sorted.bam"
#done


#===ADDED by CARO===

for c in "${CHROMOSOMES[@]}"; do
    echo "Processing $c"
    samtools view -b "$OUT_BAM" "$c" > "$CHROM_DIR/${c}.aligned.bam"
    samtools view -@ 4 -f 3 -s 42."$DECIMAL" -b "$CHROM_DIR/${c}.aligned.bam" > "$CHROM_DIR/${c}.downsampled.bam"
    samtools sort -@ 4 -o "$CHROM_DIR/${c}.downsampled.sorted.bam" "$CHROM_DIR/${c}.downsampled.bam"
    samtools index "$CHROM_DIR/${c}.downsampled.sorted.bam"
done

# === NEW STEP 4: MERGE ALL BAMS ===
samtools merge SRR23782967_downsampled.sorted.bam $CHROM_DIR/*.downsampled.sorted.bam
samtools index SRR23782967_downsampled.sorted.bam

# === OUTPUT ===
echo " Done! Output files:"
echo " - SRR23782967_downsampled.sorted.bam"
echo " - SRR23782967_downsampled.sorted.bam.bai"

# === FINAL OUTPUT ===
#echo "Done!"
#echo " - Final downsampled BAM: downsampled.sorted.bam"
#echo " - Per-chromosome BAMs: $CHROM_DIR/*.sorted.bam"

