#!/bin/bash
#SBATCH --job-name=snp_density_CAxArg
#SBATCH --output=snp_density_CAxArg_%j.out
#SBATCH --error=snp_density_CAxArg_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --partition=medmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=clazari@ucsc.edu

#====Load the bedtools module====
module load bio/bcftools
module load bio/bedtools

#==== User config ====

BCF_FILE="/share/swfsc/clazari/Omykiss/results_mega-non-model_CAxArg_March_2025/bqsr-round-0/bcf/all.bcf"     # Input .bcf file
GENOME_FILE="/share/swfsc/clazari/Omykiss/snp_density/genome_file.txt"             # Genome size file: tab-separated with 'chr<TAB>length'
WINDOW_SIZE=50000                   # Window size (e.g., 100000 for 100kb)
STEP_SIZE=10000                     # Step size (e.g., 50000 for 50kb)
OUT_PREFIX="snp_density"            # Output prefix


#----Extract SNPs from BCF file----

bcftools view -v snps "$BCF_FILE" -Oz -o snps.vcf.gz
bcftools index snps.vcf.gz


#----Create sliding window with bedtools----

bedtools makewindows -g "$GENOME_FILE" -w 50000 -s 10000 > sliding_windows_50kb_10kb.bed


#----Get SNP position in BED fromat----

bcftools query -f '%CHROM\t%POS0\t%POS\n' snps.vcf.gz > snps.bed


#----Count SNPs in each sliding window----

bedtools intersect -a sliding_windows_50kb_10kb.bed -b snps.bed -c > snp_counts_per_window.bed


#----Calculate SNP density (SNPs per bed_ per window)----

awk '{print $1, $2, $3, $4, $4/50}' snp_counts_per_window.bed > snp_density_per_50kb_window.bed


