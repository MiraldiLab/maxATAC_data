#!/bin/bash

################### generate_DHS_coverage_bigwig.sh ###################
# This workflow will window Dnase I cut site and generate a coverage track that is RP20M normalized

# INPUT:

# 1: The DHS cut sites
# example: GM12878_DHS_IS.bed.gz

# 2: Output directory
# example: ./test

# 3: Blacklisted regions
# example: hg38.composite.blacklist.bed

# 4: Slop size
# example: 10

# 5: Chromosome sizes file
# example: hg38.chrom.sizes

# 6: Scale Factor
# example: 20000000

# 7: Threads
# example: 6

# OUTPUT:

# DHS cut site bigwig track

###########################################################

### Rename Input Variables ###
DHS_bed=${1}
output_directory=${2}
blacklist=${3}
slop=${4}
chrom_sizes=${5}
scale_factor=${6}
threads=${7}

### Build names ###
base_filename=`basename ${DHS_bed} .bed.gz`
prefix=${base_filename}_RP20M_slop${slop}
bedgraph=${prefix}.bg
bigwig=${prefix}.bw

# Make directory and change into it
mkdir -p ${output_directory}

### Process ###
# http://www.metagenomics.wiki/tools/samtools/number-of-reads-in-bam-file
mapped_reads=$(zcat ${DHS_bed} | wc -l)
reads_factor=$(bc -l <<< "1/${mapped_reads}")
rpm_factor=$(bc -l <<< "${reads_factor} * ${scale_factor}")

echo "Scale factor: " ${rpm_factor}

# Use bedtools to window the signal around the cut site then remove blacklisted regions and generate coverage
# Use the total mapped reads to scale the DHS signal to RP20M, sort the bedgraph
bedtools slop -b ${slop} -g ${chrom_sizes} -i ${DHS_bed} | \
bedtools intersect -a - -b ${blacklist} -v | \
sort -k1,1 | \
bedtools genomecov -i - -g ${chrom_sizes} -bg -scale "${rpm_factor}" | \
LC_COLLATE=C sort -k1,1 -k2,2n > ${output_directory}/${bedgraph}

# Use bedGraphToBigWig to convert bedgraph to bigwig
echo "Using bedGraphToBigWig to convert bedgraph to bigwig"

bedGraphToBigWig ${output_directory}/${bedgraph} ${chrom_sizes} ${output_directory}/${bigwig}

rm "${output_directory}/${bedgraph}"