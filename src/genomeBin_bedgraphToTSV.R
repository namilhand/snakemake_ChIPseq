#!/applications/R/R-4.0.0/bin/Rscript

# This R script is called by ../Snakefile
# Convert bedgraph files into TSV files
# These TSV files can be imported into R for plotting chromosome-scale profiles 

# What does the TSV format mean in this script?
# 1. A bin that is larger than genomeBinSize (two or more successive windows with the same value are integrated into one bin in the bedgraph format) will be splitted by genomebinSize.
# 2. The last bin of the chromosome that might be smaller than genomeBinSize is filled with the value of the previous bin.

# Usage:
# Rscript genomeBin_bedgraphToTSV.R <bedg> <refgenome.fasta> <genomeBinSize>


args <- commandArgs(trailingOnly = T)
bedg_input <- args[1]
refgenome <- args[2]
genomeBinSize <- as.integer(args[3])

#bedg_input <- "results/MNase_Col_log2ChIP_binSize100Kb.bg"
#refgenome <- "/home/nison/work/refgenome/TAIR10/TAIR10.fasta"
#genomeBinSize <- 100000

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp") 
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
}

options(stringsAsFactors = F)
library(parallel)
library(plyr)
library(data.table)
library(stringr)

# Genomic definitions
fai <- read.table(paste0(refgenome, ".fai"))
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[1:5]
} else {
  chrs <- fai[,1][1:5]
}
chrLens <- fai[,2][1:5]

## Make chromosomal coordinates cumulative
## such that the first coordinate of Chr2 is
## equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

# Load bedgraph
sampleProfile <- read.table(bedg_input)
if(!grepl("Chr", fai[,1][1])) {
  sampleProfile$V1 <- paste0("Chr", sampleProfile$V1)
}
# filter out organellar genome
sampleProfile <- sampleProfile[sampleProfile$V1 %in% chrs,]

# 1. Split bins that are larger than the genomeBinSize
## Rows where the difference between end and start coordinates is > genomeBinSize
sampleProfile_bigWins <- sampleProfile[sampleProfile$V3-sampleProfile$V2 > genomeBinSize,]
## Rows where the difference between end and start coordinates is == genomeBinSize
sampleProfile <- sampleProfile[sampleProfile$V3-sampleProfile$V2 == genomeBinSize,]

## Create a list of big windows, each split into windows of genomeBinSize.
if( nrow(sampleProfile_bigWins) > 0){
    sampleProfile_bigWinsList <- mclapply(seq_along(1:dim(sampleProfile_bigWins)[1]), function(x) {
      bigWinsSplit <- seq(from = sampleProfile_bigWins[x,]$V2,
                          to = sampleProfile_bigWins[x,]$V3,
                          by = genomeBinSize)

      if(bigWinsSplit[length(bigWinsSplit)] < sampleProfile_bigWins[x,]$V3) {
        data.frame(V1 = as.character(sampleProfile_bigWins[x,]$V1),
                   V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                     bigWinsSplit[length(bigWinsSplit)])),
                   V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+genomeBinSize,
                                     sampleProfile_bigWins[x,]$V3)),
                   V4 = as.numeric(sampleProfile_bigWins[x,]$V4))
      } else if (bigWinsSplit[length(bigWinsSplit)] == sampleProfile_bigWins[x,]$V3) {
        data.frame(V1 = as.character(sampleProfile_bigWins[x,]$V1),
                   V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
                   V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+genomeBinSize),
                   V4 = as.numeric(sampleProfile_bigWins[x,]$V4))
      }
    }, mc.cores = detectCores())

    sampleProfile_bigWinsDT <- rbindlist(sampleProfile_bigWinsList)
    sampleProfile <- rbind.fill(sampleProfile, sampleProfile_bigWinsDT)
    sampleProfile <- sampleProfile[order(sampleProfile$V1, sampleProfile$V2),]
} else {
    sampleProfile <- sampleProfile[order(sampleProfile$V1, sampleProfile$V2),]
}

# 2. A bin that is smaller than genomeBinSize (at chromosome end) will be filled with the value of the previous bin.
chrLenValsList <- mclapply(seq_along(chrs), function (x) {
  chrProfileSample <- sampleProfile[sampleProfile$V1 == chrs[x],]
  if(chrProfileSample[dim(chrProfileSample)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileSample[dim(chrProfileSample)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileSample[dim(chrProfileSample)[1],]$V4))
  }
}, mc.cores = length(chrs))
sampleProfile_chrLenValsDT <- rbindlist(chrLenValsList)
sampleProfile <- rbind.fill(sampleProfile, sampleProfile_chrLenValsDT)
sampleProfile <- sampleProfile[order(sampleProfile$V1, sampleProfile$V2),]

chrCumWindowList <- lapply(seq_along(chrs), function(x) {
  chrProfileSample <- sampleProfile[sampleProfile$V1 == chrs[x],]
  chrProfileSample$V2+1+sumchr[x]
})

sampleProfile <- data.frame(chr = as.character(sampleProfile$V1),
                            window = as.integer(sampleProfile$V2+1),
                            cumwindow = as.integer(unlist(chrCumWindowList)),
                            cov = as.numeric(sampleProfile$V4))
write.table(sampleProfile,
            file = str_replace(bedg_input, "bg$", "tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)