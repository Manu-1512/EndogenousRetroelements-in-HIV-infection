#!/bin/bash
#SBATCH -p int
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 2-0:00:00
#SBATCH --mem=150G


source activate mamba_salmon

# Define the index and output directory
index="/scratch1/users/umut.cakir/transcriptome_alignment_hg38_mm10/hg38/index_files/Genes_TEindex/"
output_directory="/scratch1/users/umut.cakir/20240927_HN00226527_TRR_20240927_HN00225287_TRR/runs/Salmon/output"  # Change this to your desired output path

# Loop through all paired-end FASTQ files in the directory
for file1 in  /scratch1/users/umut.cakir/20240927_HN00226527_TRR_20240927_HN00225287_TRR/combined_data_after_adapter_removed/*_1_cleaned.fastq.gz; do
    # Generate the corresponding _2 file name
    file2="${file1/_1_cleaned.fastq.gz/_2_cleaned.fastq.gz}"

    # Check if the second file exists
    if [[ -f "$file2" ]]; then
        # Extract base name without the path and file extension
        base_name=$(basename "$file1" _1_cleaned.fastq.gz)

        # Run Salmon quantification
        salmon quant -i "$index" -l A \
        -1 "$file1" \
        -2 "$file2" \
        -o "$output_directory/$base_name" --validateMappings

        echo "Processed $file1 and $file2"
    else
        echo "Warning: Corresponding file for $file1 not found."
    fi
done

