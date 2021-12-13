#!/bin/bash

################### ATAC_pipeline.sh ###################
# Workflow to process ATAC-seq BAM. 

# Inputs:
# 1. BAM file
# 2. Sample name used for output files
# 3. Output directory
# 4. Number of cores
# 5. Blacklist .bed file
# 6. Chromosome sizes file
# 7. Flanking size or slop size to use for inferring Tn5 sites
# 8. Millions factor for normalizing signal
# 9. Deduplication flag. Whether to deduplicate or not: Choose deduplicate or skip

# Outputs:
# 1. sample_IS.bed: Insertion sites bed file
# 2. Tn5 sites file: Insertion sites that have been slopped by flanking size base pairs
# 3. RPM normalized bigwig: Reads per million normalized bigwig
#########################################################

### Inputs ###
bam=${1}
name=${2}
outDir=${3}
cores=${4}
blacklist=${5}
chromSizes=${6}
flanking_size=${7}
scale_factor=${8}
dedup=${9}

keepChr='chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22'

### Print Parameters ###
echo "Inputs:"
echo "BAM: ${bam}"
echo "Sample name: ${name}"
echo "Cores: ${cores}"
echo "Chr Keep: ${keepChr}"
echo "Blacklist: ${blacklist}"
echo "Chr Sizes: ${chromSizes}"
echo "Slop Size: ${flanking_size}"
echo "Scale factor: ${scale_factor}"

###SETUP###
deduped=${name}_deduped.bam
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
case ${dedup} in
    deduplicate)
            filtered=${name}_filtered.bam
            fix_mate=${name}_fixmate.bam

            echo "Filtering with Samtools"

            echo "Samtools sort reads by name for ${bam}"
            samtools sort -@ ${cores} -n -o ${filtered} ${bam}

            # Fixmate and sort by POSITION then index
            echo "Samtools fixmate on ${filtered}"
            samtools fixmate -@ ${cores} -m ${filtered} - | \
            samtools sort -@ ${cores} -o ${fix_mate} -
            samtools index -@ ${cores} ${fix_mate}
            rm ${filtered}

            # Mark duplicates, remove, sort, index
            echo "Remove duplicates from ${fix_mate}"
            samtools markdup -@ ${cores} -r -s ${fix_mate} - | \
            samtools sort -@ ${cores} -o ${deduped} -
            samtools index -@ ${cores} ${deduped}
            rm ${fix_mate} ${fix_mate}.bai

        ;;
    skip)
            # Sort and index file
            samtools sort -@ ${cores} -o ${deduped} ${bam}
            samtools index -@ ${cores} ${deduped}

        ;;
    *)
esac

# Remove singleton reads (-f 3) and unwanted chromosomes, sort, index. 
echo "Remove unwanted chr from ${deduped}"
samtools view -@ ${cores} -f 3 -b ${deduped} ${keepChr} | \
samtools sort -@ ${cores} -o ${final_bam} -
samtools index -@ ${cores} ${final_bam}
rm ${deduped} ${deduped}.bai

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
