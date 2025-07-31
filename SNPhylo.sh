#!/bin/bash
#SBATCH --job-name=snphylo
#SBATCH --output=snphylo.out
#SBATCH --error=snphylo.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=medmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=clazari@ucsc.edu

# -------- CONFIG --------
# Full path to your BCF input file
INPUT=/share/swfsc/clazari/Omykiss/results_mega-non-model_CAxArg_March_2025/bqsr-round-0/bcf/all.bcf

# Output prefix for the SNPhylo tree files
OUT_PREFIX=my_phylo_tree

# -------- ENV SETUP --------
conda init
conda activate snphylo_env #this has muscle and python installed
module load bio/bcftools
module load bio/samtools
module load R

# -------- RUNNING --------

# Step 1: Convert BCF to VCF in-place
bcftools view "$INPUT" -Ov -o input.vcf

# Step 2: Run SNPhylo from current directory
./snphylo.sh -v input.vcf -m 0.1 -M 0.2 -c 0.2 -o "$OUT_PREFIX"

# Done! All output (tree, SNP list, PDF, etc.) stays in this directory

