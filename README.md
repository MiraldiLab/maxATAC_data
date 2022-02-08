
![maxATAC_data_logo](https://user-images.githubusercontent.com/47329147/152259165-121fc0cd-60d2-4e6b-a72f-cbc91a095bf8.png)

# maxATAC Data

This repo contains additional resources for the maxATAC Python package.

## `./hg38`

This directory contains:

1. `hg38_maxatac_blacklist.bed`: maxATAC extended blacklist 
2. `hg38_maxatac_blacklist.bw`: 
3. `hg38.chrom.sizes`: chromosome sizes

## `./models`

Information on 127 TF models trained. Each directory is a TF name. Each directory contains a `.h5` of the best TF model, a `.tsv` file of threshold statistics, and a `.png` of the statistcs. 

Example directory structure:

```pre
|____ATF2
| |____ATF2_validationPerformance_vs_thresholdCalibration.png
| |____ATF2_binary_revcomp99_fullModel_RR0_73.h5
| |____ATF2_validationPerformance_vs_thresholdCalibration.tsv
```

## `./scripts`

A directory of scripts that are used to process the ATAC-seq, DHS, and ChIP-seq data analyzed by our publication. The scripts `ATAC_bowtie2_pipeline.sh` and `scatac_generate_bigwig.sh` are required for running `maxatac prepare`.

```pre
.
|____ATAC
| |____align_reads_STAR.sh
| |____infer_insertion_sites.sh
| |____ATAC_bowtie2_nodedup_pipeline.sh
| |____RPM_normalize_Tn5_counts.sh
| |____shift_reads.sh
| |____trim_reads.sh
| |____infer_Tn5_sites.sh
| |____scatac_generate_bigwig.sh
| |____ATAC_bowtie2_pipeline.sh
| |____sra2fastq_PE.sh
| |____filter_BAM_STAR.sh
| |____ATAC_post_alignment.sf
|____DHS
| |____infer_DHS_cut_sites.sh
| |____generate_DHS_coverage_bigwig.sh
| |____filter_DHS_SE_ENCODE.sh
|____CHIP
| |____infer_5prime_BAM.sh
| |____filter_chroms_from_bam.sh
| |____slop_5prime.sh
| |____RPM_normalize_5prime_counts.sh
```

## `Example_TF_META.tsv`

This `.tsv` file shows an example meta file required to train your own models. 

## Publication

The maxATAC pre-print is currently available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.01.28.478235v1.article-metrics). 

```pre
maxATAC: genome-scale transcription-factor binding prediction from ATAC-seq with deep neural networks
Tareian Cazares, Faiz W. Rizvi, Balaji Iyer, Xiaoting Chen, Michael Kotliar, Joseph A. Wayman, Anthony Bejjani, Omer Donmez, Benjamin Wronowski, Sreeja Parameswaran, Leah C. Kottyan, Artem Barski, Matthew T. Weirauch, VB Surya Prasath, Emily R. Miraldi
bioRxiv 2022.01.28.478235; doi: https://doi.org/10.1101/2022.01.28.478235
```
