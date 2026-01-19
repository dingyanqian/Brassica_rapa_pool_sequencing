#!/bin/bash
# This script annotates genomic regions for selection analysis
# For each treatment, it performs the following steps:
# 1. Creates proper BED format for regions
# 2. Uses bedtools intersect to annotate regions with the reference GFF
# 3. Extracts gene names from the annotated regions using awk
# 4. Links significant genes to their functional annotations using grep

# Input file with genomic regions
input_file="../1_input/1_selected_regions_summary.txt"

# Reference files (adjust paths if necessary)
reference_gff="../1_input/1_reference_fixed.gff"
annotation_info="../1_input/1_BrapaFPsc_277_v1.3.annotation_info.txt"

# Get unique treatment (skip header if present)
treatments=$(tail -n +2 "$input_file" | awk '{print $5}' | sort | uniq)

# Loop through each treatment
for treatment in $treatments; do
  echo "Processing $treatment..."
  
  # Create proper BED format (chrom, start, end) for bedtools
  awk -v t="$treatment" -v OFS='\t' '$5 == t {print $1, $2, $3}' "$input_file" > "${treatment}_regions.bed"
  
  # Step 1: Annotate regions with bedtools intersect
  bedtools intersect -a "${treatment}_regions.bed" -b "$reference_gff" -wa -wb > "${treatment}_annotated_regions.txt"
  
  # Step 2: Extract gene names using grep and awk
  # First find lines with "gene" feature type
  grep "gene" "${treatment}_annotated_regions.txt" | 
    # Then extract gene names using awk
    awk -F'\t' '{
      for(i=1; i<=NF; i++) {
        if($i ~ /Name=/) {
          split($i, a, "Name=");
          split(a[2], b, ";");
          print b[1];
          break;
        }
      }
    }' | sort | uniq > "${treatment}_annotated_genes.txt"
  
  # Step 3: Link significant genes to their functional annotations using grep
  grep -wFf "${treatment}_annotated_genes.txt" "$annotation_info" > "${treatment}_annotated_significant_genes.txt"
  
  echo "Finished processing $treatment."
  echo "  - Regions BED file: ${treatment}_regions.bed"
  echo "  - Annotated Regions: ${treatment}_annotated_regions.txt"
  echo "  - Annotated Genes: ${treatment}_annotated_genes.txt"
  echo "  - Significant Gene Annotations: ${treatment}_annotated_significant_genes.txt"
done

echo "All categories processed."