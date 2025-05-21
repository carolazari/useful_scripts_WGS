#!/bin/bash
#SBATCH -p medmem
#SBATCH --job-name=downsample_genbank_seq_parallel
#SBATCH --output=downsample_genbank_seq.log
#SBATCH --error=run-single-error-%j
#SBATCH --time=160:00:00
#SBATCH -c 10
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=clazari@ucsc.edu
# Load modules
module load bio/samtools
module load aligners/bwa

set -e

# === INPUT ===
#R1="/home/clazari/Projects/Chinook/mega-non-model-wgs-snakeflow/fastq_run2/fastq/SRR23782967_1.fastq"
#R2="/home/clazari/Projects/Chinook/mega-non-model-wgs-snakeflow/fastq_run2/fastq/SRR23782967_2.fastq"
R1="../fastq_run2/fastq/SRR23782967_1.fastq"
R2="../fastq_run2/fastq/SRR23782967_2.fastq"
REF="../resources/genome.fasta"
#REF="~/refs/Omykiss_Arlee_genome/genome.fasta"
GENOME_SIZE=2406167875   # e.g. 1 Gb (adjust to your genome)
READ_LENGTH=150          # e.g. 150 bp
TARGET_DEPTH=8

# === PREP ===
# Ensure reference is indexed
bwa index $REF
samtools faidx $REF

# === STEP 1: ALIGN PAIRED-END READS ===old version
#bwa mem -t 4 $REF $R1 $R2 | samtools view -b -o aligned.bam -

# === STEP 1: ALIGN PAIRED-END READS ===
echo "[INFO] Aligning paired-end reads..."
bwa mem -t 4 $REF $R1 $R2 | samtools view -b -o aligned.bam -

# === CHECK BAM INTEGRITY ===
echo "[INFO] Checking BAM file integrity..."
if ! samtools quickcheck -v aligned.bam; then
    echo "[WARNING] aligned.bam is incomplete or corrupted. Regenerating..."
    rm -f aligned.bam
    bwa mem -t 4 $REF $R1 $R2 | samtools view -b -o aligned.bam -
    echo "[INFO] Rechecking aligned.bam..."
    if ! samtools quickcheck -v aligned.bam; then
        echo "[ERROR] aligned.bam is still invalid after regeneration. Exiting."
        exit 1
    fi
    echo "[INFO] aligned.bam was successfully regenerated and passed integrity check."
else
    echo "[INFO] aligned.bam passed integrity check."
fi

# === STEP 2: CALCULATE CURRENT COVERAGE ===
TOTAL_READS=$(samtools view -c aligned.bam)
CURRENT_DEPTH=$(echo "$TOTAL_READS * $READ_LENGTH / $GENOME_SIZE" | bc -l)
FRACTION=$(echo "$TARGET_DEPTH / $CURRENT_DEPTH" | bc -l)

echo "Total reads: $TOTAL_READS"
echo "Estimated depth: $CURRENT_DEPTH"
echo "Sampling fraction for ~${TARGET_DEPTH}×: $FRACTION"

# === STEP 3: DOWNSAMPLE, SORT, AND INDEX BAM ===
# Extract only the decimal part for samtools -s argument
DECIMAL=$(printf "%.3f" $FRACTION | cut -d'.' -f2)

#=== OLD VERSION ===
#VAR=$( samtools idxstats aligned.bam | cut -f 1 )
#readarray -t CHROMOSOMES <<<"$VAR"

#===ADDED by CARO===
# Get chromosome names into a Bash array
readarray -t CHROMOSOMES < <(samtools idxstats aligned.bam | cut -f1 | grep -v '\*')

declare -p CHROMOSOMES

#=== OLD ===

#for[ c in $CHROMOSOMES ]; do
#    samtools view -b aligned.bam $c > ${DIR}/ ${c}.aligned.bam
#    samtools view -@ 4 -f 3 -s 42.$DECIMAL -b ${DIR}/${c}.aligned.bam > ${DIR}/${c}.downsampled.bam
#    samtools sort -@ 4 -o ${DIR}/${c}.downsampled.sorted.bam ${DIR}/${c}.downsampled.bam
#    samtools index ${DIR}/${c}.downsampled.sorted.bam
#done

#===ADDED by CARO===

for c in "${CHROMOSOMES[@]}"; do
    echo "Processing $c"
    samtools view -b aligned.bam "$c" > "${DIR}/${c}.aligned.bam"
    samtools view -@ 4 -f 3 -s 42."$DECIMAL" -b "${DIR}/${c}.aligned.bam" > "${DIR}/${c}.downsampled.bam"
    samtools sort -@ 4 -o "${DIR}/${c}.downsampled.sorted.bam" "${DIR}/${c}.downsampled.bam"
    samtools index "${DIR}/${c}.downsampled.sorted.bam"
done

# === NEW STEP 4: MERGE ALL BAMS ===
samtools merge SRR23782967_downsampled.sorted.bam ${DIR}/*.downsampled.sorted.bam
samtools index SRR23782967_downsampled.sorted.bam

# === OUTPUT ===
echo " Done! Output files:"
echo " - SRR23782967_downsampled.sorted.bam"
echo " - SRR23782967_downsampled.sorted.bam.bai"

###  ==== OLD CODE ===
### samtools view -@ 4 -f 3 -s 42.$DECIMAL -b aligned.bam > downsampled.bam

# === STEP 4: SORT AND INDEX ===
## samtools sort -@ 4 -o SRR23782967_downsampled.sorted.bam SRR23782967_downsampled.bam
## samtools index SRR23782967_downsampled.sorted.bam