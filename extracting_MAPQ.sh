#!/bin/bash
set -euo pipefail

mkdir -p mapq_hists
echo "sample,mapq,count" > combined_mapq_counts.csv

for bam in *.bam; do
  sample="${bam%.bam}"
  echo "Processing $sample (chromosome NC_048571.1)..."

  samtools view "$bam" NC_048571.1 \
    | awk '{print $5}' \
    | sort \
    | uniq -c \
    | awk -v s="$sample" '{print s "," $2 "," $1}' \
    > mapq_hists/"$sample"_NC_048571.1_mapq.csv

  cat mapq_hists/"$sample"_NC_048571.1_mapq.csv >> combined_mapq_counts.csv
done

echo "Done! Individual MAPQ CSVs in mapq_hists/, combined CSV in combined_mapq_counts.csv"

