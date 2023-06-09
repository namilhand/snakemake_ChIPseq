# Average each bins of deeptools:computeMatrix tab output, then covert to tsv format
# input: deeptools:computeMatrix tab output

library(tidyverse)

args <- commandArgs(trailingOnly=T)
input.file <- args[1]
upstream <- as.numeric(args[2])
downstream <- as.numeric(args[3])
binSize <- as.numeric(args[4])
output <- args[5]

#input.file <- "/datasets/data_4/nison/NGS_libraries/MNase/Choi18_WT_suf3-1/results/05_metaprofile/genes/MNase_suf3-bud_metaprofile_region-genes_referencePoint-TES.tab"
#upstream <- 2000
#downstream <- 2000
#binSize <- 10
#output <- "test.tsv"

nbin <- (upstream + downstream)/binSize
binName <- paste0("bin", 1:nbin)

input <- read_tsv(input.file, skip=3, col_names=F)
colnames(input) <- binName

# Calculate average score of each bin
avgScore.tsv <- input %>%
    summarise(across(.cols=1:ncol(input), ~ mean(.x, na.rm=TRUE))) %>%
    pivot_longer(cols = 1:ncol(input), names_to="bin", values_to="meanScore")

write_tsv(avgScore.tsv, file=output, col_names=T)
