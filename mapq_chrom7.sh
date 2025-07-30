#!/bin/bash
#SBATCH --job-name=mapq_summary_caxnz
#SBATCH --output=mapq_summary_%j.out
#SBATCH --error=mapq_summary_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=medmem
#SBATCH --mail-user=clazari@ucsc.edu

# Load modules (update to match your cluster's setup)
module load bio/samtools
module load bio/bedtools

# Input and output setup
GENOME_FA=/home/clazari/Projects/Omykiss/mega-non-model-wgs-snakeflow/resources/genome.fasta               # <-- path to your reference genome
#BAM_DIR=/share/swfsc/clazari/Omykiss/results_mega-non-model_CAxArg_March_2025/bqsr-round-0/overlap_clipped/ 
BAM_DIR=/share/swfsc/clazari/Omykiss/results_mega-non-model_CAxNZ_March25/bqsr-round-0/overlap_clipped/               # <-- path to your .bam files
OUT_DIR=mapq_genome-wide        # <-- output directory

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Generate genome index and windows (if not already done)
if [ ! -f "$GENOME_FA.fai" ]; then
    samtools faidx "$GENOME_FA"
fi

bedtools makewindows -g "$GENOME_FA.fai" -w 50000 > genome_windows.bed

# Loop over each BAM and process
for bam in "$BAM_DIR"/*.bam; do
    sample=$(basename "$bam" .bam)
    echo "Processing $sample"

    # Extract MAPQ data as BED format
    samtools view "$bam" | awk -v OFS='\t' '{print $3, $4-1, $4, $5}' > "${sample}_mapq.bed"

    # Map MAPQ values to genome windows (average per window)
    bedtools map -a genome_windows.bed -b "${sample}_mapq.bed" -c 4 -o mean \
        > "$OUT_DIR/${sample}_mapq_windows.bed"

    # Clean up intermediate
    rm "${sample}_mapq.bed"
done
