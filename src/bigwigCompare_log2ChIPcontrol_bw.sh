#!/bin/bash

# Usage:
# csmit -m 50G -c 48 "bash ./bigwigCompare_log2ChIPcontrol_bw.sh WT_SPO11oligos_Rep1 WT_gDNA_Rep1_R1 1000 1kb 48"

ChIPName=$1
controlName=$2
binSize=$3
binName=$4
threads=$5

[ -d ../mapped/both/log2ChIPcontrol ] || mkdir -p ../mapped/both/log2ChIPcontrol

source activate ChIPseq_mapping

bigwigCompare -b1 "../mapped/both/bw/"${ChIPName}"_MappedOn_TAIR10_chr_all_lowXM_both_sort_norm.bw" \
              -b2 "/home/ajt200/analysis/150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos_TAIR10_chr_all/mapped/both/bw/"${controlName}"_MappedOn_TAIR10_chr_all_lowXM_both_sort_norm.bw" \
              -o "../mapped/both/log2ChIPcontrol/log2_"${ChIPName}"_"${controlName}"_MappedOn_TAIR10_chr_all_lowXM_both_sort_norm_binSize"${binName}".bw" \
              -of "bigwig" \
              --pseudocount 1 \
              --operation log2 \
              --binSize ${binSize} \
              -p ${threads}

source deactivate
