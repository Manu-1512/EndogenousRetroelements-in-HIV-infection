# ====================================
# Author:Manu Singh 20 November 2020

# =====================================
# Set Enviornment variables as per the .bash_profiles on the system
# =====================================

#first we are dealing with TE coordinates of human genome
# the file name is  "Human_class2_TE.bed"
#Human_class2_TE.bed file looks like following

#chr1	41379	42285	ERVL-E-int	1118.000000	-	hg19_rmsk	exon	.	gene_id "ERVL-E-int"; transcript_id "ERVL-E-int_dup3"; 

# these bed files are already curated based on the 3 step filtering process mentioned in manuscript
# this file was generated and further modified for the analysis of manuscripts viz.  Nature 516 (7531), 405-409, Bioessays 38 (1), 109-117, Circulation 136 (19), 1824-1839, JCI insight 5 (7), Viruses 12 (11), 1303, bioRxiv, 318329
# the full bed files are available on request 

#Here we are going to  make a combined gtf file containing these TEs and gene from human hg19 sequences

# First step to create an ID for each sequences, which we could fetch later for the integrative kinds of analysis
 
cut -f 2 Human_class2_TE.bed -d ";" | cut -f3 -d " " | sed 's/"//g' > Human_class2_TE_transID.txt

## Here, we add the prefix "TELONG", so that we can fetch TEs from mixed data later

paste Human_class2_TE.bed Human_class2_TE_transID.txt | awk '{print $1,$2,$3,"TELONG-"$14,$5,$6}' OFS="\t" > Human_class2_TE_ID.bed

awk '{print $1,$2,$3,$4"-"$1"-"$2"-"$3,$5,$6}' OFS="\t" Human_class2_TE_ID.bed > Human_class2_TE_IDv2.bed

bedtools getfasta -fi hg19.fa -bed Human_class2_TE_IDv2.bed -s -name | awk '{ print toupper($0) }' | sed 's/(+)\|(-)//g'  > Human_class2_TE_ID.fa

bedtools maskfasta -fi hg19.fa -bed Human_class2_TE_ID.bed -fo hg19_masked.fa


#above commands have given the fasta files of TEs , where the headers can be split to gain any information from TEs as given above



# Now let's deal with genes 
#This is mart export file, downlaoded from UCSC, with following Genes coordinates
#file name is "mart_export.txt" 
#ENSG00000229483	13	23743974	23744736	-1	LINC00362


awk 'NR>1{if($5=="1") print "chr"$2,$3,$4,$6"_"$1,$1,"+"; else print "chr"$2,$3,$4,$6"_"$1,$1,"-"}' OFS="\t" mart_export.txt  > hg19_Genes.bed



bedtools getfasta -fi hg19_masked.fa -bed hg19_Genes.bed -s -name | awk '{ print toupper($0) }' | sed 's/(+)\|(-)//g'  > hg19_Genes.fa

cat hg19_Genes.fa Human_class2_TE_ID.fa > hg19_Genes_long_TE_ID.fa

grep ">" hg19_Genes_long_TE_ID.fa | sed 's/>//g' | awk '{print $1,$1}' OFS="\t" > txp2gene.tsv


/programs/salmon-1.2.1/bin/salmon index -t hg19_Genes_long_TE_ID.fa -i Genes_TEindex -k 31

########