#!/bin/bash
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 2-0:00:00
#SBATCH --mem=250G
 
mkdir Repeat_analysis
cp *ChIP_summits.bed Repeat_analysis/
cd Repeat_analysis
#Repeats.bed gencode.vM10.annotation.gtf GRCm38.primary_assembly.genome.fa mart_export_genes.tsv 
wget -O "Repeats_file.bed" "https://usegalaxy.eu/api/datasets/4838ba20a6d86765aba88a3937c9f502/display?to_ext=bed" #Repeats.bed #Mammal, Human, Dec. 2013 GRCh38, Repeats, RepeatMasker, rmsk

sort -k1,1 -k2,2n Repeats_file.bed > Repeats.bed

awk '$5 != "Simple_repeat" && $5 != "Low_complexity" && $5 != "DNA" && $5 != "Satellite"' Repeats.bed | sort -k1,1 -k2,2n > Repeats_f1.bed
awk '$6 != "Alu" ' Repeats_f1.bed | sort -k1,1 -k2,2n > Repeats_f2.bed
awk '($3-$2) >=100' Repeats_f2.bed | sort -k1,1 -k2,2n > Repeats_f3.bed #any TE sequence that is < 100 bps


# Input files
repeats="Repeats_f3.bed"
bed_files=($(find . -maxdepth 1 -name '*ChIP_summits.bed' | sort))

# Temporary output files
non_binary_output="non_binary_overlap_counts.tmp"
binary_output="binary_overlap_counts.tmp"

# Copy the original repeats file as a starting point
cp "$repeats" "$non_binary_output"
cp "$repeats" "$binary_output"

# Initialize the header line with the first 6 column names
header_line="#genoName\tgenoStart\tgenoEnd\trepName\trepClass\trepFamily"

# Process each BED file for overlaps
for bedf in "${bed_files[@]}"; do
    echo $bedf
    bedtools intersect -c -a "$repeats" -b "$bedf" > temp_overlap.txt

    # Add the new non-binary counts as a new column (starting from column 7)
    paste "$non_binary_output" <(cut -f7 temp_overlap.txt) > temp_combined.txt
    mv temp_combined.txt "$non_binary_output"

    # Convert counts to binary and add as a new column
    awk '{if ($NF > 0) print 1; else print 0}' temp_overlap.txt > temp_binary_col.txt
    paste "$binary_output" temp_binary_col.txt > temp_combined.txt
    mv temp_combined.txt "$binary_output"

    # Add the current bed_file name to the header line inside the loop
    bn=$(basename "$bedf")
    header_line="${header_line}\t${bn}"
done

# Add the Total column at the end of the header line
header_line="${header_line}\tTotal"

# Add total overlaps (non-binary): sum from column 7 to the end
awk 'BEGIN{OFS="\t"}{sum=0; for(i=7;i<=NF;i++) sum+=$i; print $0, sum}' "$non_binary_output" \
    | sort -k1,1 -k2,2n > final_non_binary_output.bed

# Add total overlaps (binary): sum from column 7 to the end
awk 'BEGIN{OFS="\t"}{sum=0; for(i=7;i<=NF;i++) sum+=$i; print $0, sum}' "$binary_output" \
    | sort -k1,1 -k2,2n > final_binary_output.bed

# Cleanup
rm temp_overlap.txt temp_binary_col.txt "$non_binary_output" "$binary_output"

# Prepend the completed header line to the final output files
{ echo -e "$header_line"; cat final_non_binary_output.bed; } > temp && mv temp final_non_binary_output.bed
{ echo -e "$header_line"; cat final_binary_output.bed; } > temp && mv temp final_binary_output.bed

echo "final_non_binary_output.bed and final_binary_output.bed have been created."

wget -O "gencode.v46.annotation.gtf.gz" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz" #gencode.vM10.annotation.gtf
gunzip "gencode.v46.annotation.gtf.gz"
awk '$3 == "exon" || $3 == "UTR" {print $1 "\t" $4-1 "\t" $5 "\t" $3 "\t" $10}' gencode.v46.annotation.gtf | sed 's/;.*//g' | sed 's/"//g' | sort -k1,1 -k2,2n > exons_and_UTRs.bed
awk '$3 == "gene" {print $1 "\t" $4-1 "\t" $5 "\t" $3 "\t" $10$14}' gencode.v46.annotation.gtf | sed 's/";"/_/g; s/"/ /g; s/;/_/g; s/_$//' | sort -k1,1 -k2,2n > genes.bed

bedtools closest -a final_binary_output.bed -b genes.bed -header -d | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"} NR==1 {
   # Store original header columns in h
   h=$1
   for(i=2;i<=NF;i++) h=h OFS $i
   # Print original header plus the 6 new column names
   print h,"Closest_chromosome","Closest_genoStart","Closest_genoEnd","Closest_type","Closest_GeneID_GeneName","Closest_Distance"
   next
}1'  > final_binary_output_with_closest_genes.bed

bedtools closest -a final_non_binary_output.bed -b genes.bed -header -d | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"} NR==1 {
   # Store original header columns in h
   h=$1
   for(i=2;i<=NF;i++) h=h OFS $i
   # Print original header plus the 6 new column names
   print h,"Closest_chromosome","Closest_genoStart","Closest_genoEnd","Closest_type","Closest_GeneID_GeneName","Closest_Distance"
   next
}1' > final_non_binary_output_with_closest_genes.bed

bedtools subtract -a final_binary_output_with_closest_genes.bed -b exons_and_UTRs.bed -A -header | sort -k1,1 -k2,2n > final_binary_output_with_closest_genes_overlapping_with_exons_removed.bed #Filtered out overlapping with exons and UTRs of genes
bedtools subtract -a final_non_binary_output_with_closest_genes.bed -b exons_and_UTRs.bed -A -header | sort -k1,1 -k2,2n > final_non_binary_output_with_closest_genes_overlapping_with_exons_removed.bed #Filtered out overlapping with exons and UTRs of genes

awk -v c=$(head -1 final_binary_output_with_closest_genes_overlapping_with_exons_removed.bed | tr '\t' '\n' | nl | awk '$2=="Total"{print $1}') 'NR==1 || (NR>1 && $c>=4)' final_binary_output_with_closest_genes_overlapping_with_exons_removed.bed > final_binary_output_with_closest_genes_overlapping_with_exons_removed_detected_at_least_4_samples.bed

awk -v c=$(head -1 final_binary_output_with_closest_genes.bed | tr '\t' '\n' | nl | awk '$2=="Total"{print $1}') 'NR==1 || (NR>1 && $c>=1)' final_binary_output_with_closest_genes.bed > final_binary_output_with_closest_genes_detected_at_least_1_sample.bed

awk -v c=$(head -1 final_binary_output_with_closest_genes_overlapping_with_exons_removed.bed | tr '\t' '\n' | nl | awk '$2=="Total"{print $1}') 'NR==1 || (NR>1 && $c>=1)' final_binary_output_with_closest_genes_overlapping_with_exons_removed.bed > final_binary_output_with_closest_genes_overlapping_with_exons_removed_at_least_1_sample.bed

echo "Completed"
 
 