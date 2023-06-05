#!/bin/bash

# input bw
bw_test="/datasets/data_4/nison/NGS_libraries/ChIP_K9me2/results/03_bamCoverage/bw/ERR3813867_MappedOn_tair10_nuclear_sort_md_norm_1bp-resolution.bw"
bw_ctrl="/datasets/data_4/nison/NGS_libraries/ChIP_K9me2/results/03_bamCoverage/bw/ERR3813868_MappedOn_tair10_nuclear_sort_md_norm_1bp-resolution.bw"

# output directory
dirout="/datasets/data_4/nison/NGS_libraries/ChIP_K9me2/results/04_log2ChIP"


# misc
threads=30
toTSV="genomeBin_bedgraphToTSV.R"
refgenome_fasta="/home/nison/work/refgenome/TAIR10/TAIR10.fasta"
#======================#

mkdir -p $dirout

# turn off toTSV module in the bwcompare function if binSize = 1bp
function bwcompare {
    chip=$1
    control=$2
    prefix=$3
	binSize=$4
	binName=$5
	dirout=$6
	# bedg file will be transformed to tsv format by default. Set $7 as 0 if not desired.
	toTSV=${7:-1}

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
	if [[ $toTSV ]]; then
		Rscript $toTSV $dirout/${prefix}_log2ChIP_binSize${binName}.bg $refgenome_fasta ${binSize}
	fi
}

bwcompare $bw_test $bw_ctrl K9me2_WT-bud 1 1bp $dirout

