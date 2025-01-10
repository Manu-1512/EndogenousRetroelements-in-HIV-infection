#!/bin/bash
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 1-00:00:00
#SBATCH --mem=100G


for i in *.fastq.gz; \
      do \
      echo "$i"
      bowtie2 \
                --phred33 -p 24 \
                --very-sensitive-local \
                -x human_index/hg38 \
                -U "$i" \
                -S ${i%.fastq.gz}.sam; \
      done
 
echo "Completed"
 
