#!/bin/bash
#SBATCH --job-name=blast_regions
#SBATCH --output=blast_regions_%j.out
#SBATCH --error=sblast_regions_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=blast_regions_%j.log
#SBATCH --partition=medmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=clazari@ucsc.edu

# --------------------------------------
# SLURM job: Run BLAST for 3 genomic regions
# against the reference genome
# --------------------------------------

# Load BLAST+ module
module load bio/blast

# Input genome FASTA file
GENOME_FA="/home/clazari/Projects/Omykiss/mega-non-model-wgs-snakeflow/resources/genome.fasta"
DB_NAME="my_genome_db"

# Check if BLAST database exists, if not, create it
if [[ ! -f ${DB_NAME}.nin || ! -f ${DB_NAME}.nsq ]]; then
  echo "Creating BLAST database from $GENOME_FA..."
  makeblastdb -in "$GENOME_FA" -dbtype nucl -out "$DB_NAME"
  echo "BLAST database created."
else
  echo "BLAST database $DB_NAME already exists. Skipping makeblastdb."
fi

# Define query FASTA files
QUERIES=("NC_048571_34-34.5Mb.fa" "NC_048571_35-35.5Mb.fa" "NC_048571_36-36.5Mb.fa")

# Loop through each query and run blastn
for QUERY in "${QUERIES[@]}"
do
  BASENAME=$(basename "$QUERY" .fa)
  OUTPUT="${BASENAME}_blast.out"

  echo "Running BLAST for $QUERY..."
  blastn -query "$QUERY" \
         -db "$DB_NAME" \
         -out "$OUTPUT" \
         -evalue 1e-10 \
         -outfmt 6 \
         -max_target_seqs 10 \
         -num_threads 8

  echo "Finished BLAST for $QUERY → $OUTPUT"
done

echo "All BLAST jobs completed successfully."
