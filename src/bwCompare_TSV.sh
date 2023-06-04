#!/bin/bash

# 1. input test bw files
dir_bw_chip="/datasets/data_4/nison/NGS_libraries/MNase/Choi18_WT_suf3-1/results/03_bamCoverage/bw"
bw_col="$dir_bw_chip/mnase-col0_MappedOn_tair10_nuclear_sort_md_norm_1bp-resolution.bw"
bw_arp6="$dir_bw_chip/mnase-arp6_MappedOn_tair10_nuclear_sort_md_norm_1bp-resolution.bw"

# 2. control bw files
dir_bw_control="/datasets/data_4/nison/NGS_libraries/MNase/Choi18_MNase-control/results/03_bamCoverage/bw"
bw_control="$dir_bw_control/ERR2215860_MappedOn_tair10_nuclear_sort_md_norm_1bp-resolution.bw"

# 3. genomeBinSize
binSize=10000
binName="10Kb"

# smoothLength=300000
# smoothName="300Kb"

# 4. output directories
dirout="/datasets/data_4/nison/NGS_libraries/MNase/Choi18_WT_suf3-1/results/04_log2ChIP/results"
dirlog="/datasets/data_4/nison/NGS_libraries/MNase/Choi18_WT_suf3-1/results/04_log2ChIP/logs"

# 5. misc
threads=30
toTSV="genomeBin_bedgraphToTSV.R"
refgenome_fasta="/home/nison/work/refgenome/TAIR10/TAIR10.fasta"
#======================#

mkdir -p $dirout
mkdir -p $dirlog

function bwcompare {
    chip=$1
    control=$2
    prefix=$3

	log2_output=${prefix}_log2ChIP_binSize${binName}.${suffix}


    bigwigCompare -b1 $chip -b2 $control -of bedgraph \
        --binSize $binSize \
        -p $threads \
		--pseudocount 1 \
		--operation log2 \
        --skipZeroOverZero \
        -o $dirout/${prefix}_log2ChIP_binSize${binName}.bg
    # --skipZeroOverZero: Skip bins where BOTH BAM files lack coverage.
    bigwigCompare -b1 $chip -b2 $control -of bigwig \
        --binSize $binSize \
        -p $threads \
		--pseudocount 1 \
		--operation log2 \
        --skipZeroOverZero \
        -o $dirout/${prefix}_log2ChIP_binSize${binName}.bw

	Rscript $toTSV $dirout/${prefix}_log2ChIP_binSize${binName}.bg $refgenome_fasta ${binSize}
}

bwcompare $bw_col $bw_control MNase_col
bwcompare $bw_col $bw_control MNase_arp6

