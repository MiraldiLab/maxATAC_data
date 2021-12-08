#!/bin/bash

################### ATAC_pipeline.sh ###################
# Workflow to process ATAC-seq fastqs. This pipeline uses bowtie2 to map reads

# Inputs:
# 1. BAM file
# 2. Sample name used for output files
# 3. Output directory
# 4. Number of cores
# 5. Blacklist .bed file
# 6. Chromosome sizes file
# 7. Flanking size or slop size to use for inferring Tn5 sites
# 8. Millions factor for normalizing signal

#########################################################

###Inputs###
bam=${1}
name=${2}
outDir=${3}
cores=${4}
keepChr='chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22'
blacklist=${5}
chromSizes=${6}
flanking_size=${7}
scale_factor=${8}

echo "Inputs:"
echo "BAM: ${bam}"
echo "Sample name: ${name}"
echo "Cores: ${cores}"
echo "Bowtie2 index: ${idxBowtie2}"
echo "Chr Keep: ${keepChr}"
echo "Blacklist: ${blacklist}"
echo "Chr Sizes: ${chromSizes}"
echo "Slop Size: ${flanking_size}"
echo "Millions factor: ${millions_factor}"

###SETUP###
final_bam=${name}_final.bam
insertion_sites=${name}_IS.bed.gz
Tn5_sites=${name}_IS_slop${flanking_size}.bed.gz
bedgraph=${name}_IS_slop${flanking_size}_RP20M.bg
bigwig=${name}_IS_slop${flanking_size}_RP20M.bw

# make and move to the working directory
echo "Create working directory: ${outDir}"
mkdir -p ${outDir}
cd ${outDir}

###Process###
echo "Filtering with Samtools"
# Remove singleton reads (-f 3) and unwanted chromosomes, sort, index. 
echo "Remove unwanted chr from ${bam}"
samtools view -@ ${cores} -f 3 -b ${bam} ${keepChr} | \
samtools sort -@ ${cores} -o ${final_bam} -
samtools index -@ ${cores} ${final_bam}

# Infer insertion sites at 1 bp resolution
bedtools bamtobed -i ${final_bam} | \
awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $2 + 5, $4, $5, $6; else print $1, $3 - 5, $3 - 4, $4, $5, $6}' | \
bedtools intersect -a - -b ${blacklist} -v | \
pigz > ${insertion_sites}

# Bedtools can take as input a compressed or uncompressed BED file. 
bedtools slop  -i ${insertion_sites} -g ${chromSizes} -b ${flanking_size} | \
bedtools intersect -a - -b ${blacklist} -v | \
pigz > ${Tn5_sites}

echo "Scale factor: " "${scale_factor}"

# Use bedtools to obtain genome coverage
echo "Using Bedtools to convert BAM to bedgraph"
bedtools genomecov -i ${Tn5_sites} -g ${chromSizes} -bg -scale ${scale_factor} | LC_COLLATE=C sort -k1,1 -k2,2n > ${bedgraph}

# Use bedGraphToBigWig to convert bedgraph to bigwig
echo "Using bedGraphToBigWig to convert bedgraph to bigwig"

bedGraphToBigWig ${bedgraph} ${chromSizes} ${bigwig}

rm ${bedgraph}

echo 'Done!'
