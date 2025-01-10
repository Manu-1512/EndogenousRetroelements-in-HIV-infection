#!/bin/bash
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 0-12:00:00
#SBATCH --mem=90G



wget "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

mkdir hg38
mv hg38.fa.gz hg38/
cd hg38/
gzip -d hg38.fa.gz

cd ..
mkdir human_index

bowtie2-build hg38/hg38.fa human_index/hg38