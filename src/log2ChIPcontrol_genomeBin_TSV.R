#!/applications/R/R-4.0.0/bin/Rscript

# Calculate unsmoothed and moving-average-smoothed
# log2(ChIP/control) chromosome-scale profiles 

# Usage:
# ./log2ChIPcontrol_genomeBin_TSV.R WT_SPO11oligos_Rep1 WT_gDNA_Rep1_R1 TAIR10_chr_all both 10000 101

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
controlName <- args[2]
refbase <- args[3]
align <- args[4]
genomeBinSize <- as.integer(args[5])
maPeriod <- as.integer(args[6])

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

outDir <- "log2ChIPcontrol/"
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))

# Genomic definitions
fai <- read.table(paste0("../../../data/index/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1][1:5])
} else {
  chrs <- fai[,1][1:5]
}
chrLens <- fai[,2][1:5]

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

# Load bedgraph
sampleProfile <- read.table(paste0(sampleName, "_MappedOn_", refbase, "_lowXM_",
                                   align, "_sort_norm_binSize", genomeBinName, ".tsv"),
                            header = T)
if(controlName == "WT_gDNA_Rep1_R1") {
  controlProfile <- read.table(paste0("/home/ajt200/analysis/150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos_", refbase, "/mapped/",
                                      align, "/tsv/",
                                      controlName, "_MappedOn_", refbase, "_lowXM_",
                                      align, "_sort_norm_binSize", genomeBinName, ".tsv"),
                               header = T)
} else {
  if(!(controlName %in% c("WT_gDNA_Rep1_R1"))) {
    stop("controlName is not WT_gDNA_Rep1_R1")
  }
}

covProfile <- data.frame(chr = sampleProfile$chr,
                         window = sampleProfile$window,
                         cumwindow = sampleProfile$cumwindow,
                         ChIP = sampleProfile[,4],
                         control = controlProfile[,4],
                         log2ChIPcontrol = log2( (sampleProfile[,4] + 1) /
                                                 (controlProfile[,4] + 1) ))

chr_covProfile <- lapply(seq_along(chrs), function(x) {
  covProfile[covProfile$chr == chrs[x],]
})

# Calculate moving average of current window,
#### (maPeriod/2) previous windows (where maPeriod is even) OR
# (maPeriod/2)-0.5 previous windows (where maPeriod is odd),
# and
#### (maPeriod/2) subsequent windows (where maPeriod is even) OR
# (maPeriod/2)-0.5 subsequent windows (where maPeriod is odd)
# (the higher maPeriod is, the greater the smoothing)
stopifnot(maPeriod %% 2 != 0)
flank <- (maPeriod/2)-0.5
# Define MA filter coefficients
f <- rep(1/maPeriod, maPeriod)

filt_chr_covProfile <- mclapply(seq_along(chr_covProfile), function(x) {
  # log2ChIPcontrol
  filt_log2ChIPcontrol <- stats::filter(x = chr_covProfile[[x]]$log2ChIPcontrol,
                                        filter = f,
                                        sides = 2)
  filt_log2ChIPcontrol[1:flank] <- filt_log2ChIPcontrol[flank+1]
  filt_log2ChIPcontrol[(length(filt_log2ChIPcontrol)-flank+1):length(filt_log2ChIPcontrol)] <- filt_log2ChIPcontrol[(length(filt_log2ChIPcontrol)-flank)]

  # ChIP
  filt_ChIP <- stats::filter(x = chr_covProfile[[x]]$ChIP,
                             filter = f,
                             sides = 2)
  filt_ChIP[1:flank] <- filt_ChIP[flank+1]
  filt_ChIP[(length(filt_ChIP)-flank+1):length(filt_ChIP)] <- filt_ChIP[(length(filt_ChIP)-flank)]

  # control
  filt_control <- stats::filter(x = chr_covProfile[[x]]$control,
                                filter = f,
                                sides = 2)
  filt_control[1:flank] <- filt_control[flank+1]
  filt_control[(length(filt_control)-flank+1):length(filt_control)] <- filt_control[(length(filt_control)-flank)]

  data.frame(chr = as.character(chr_covProfile[[x]]$chr),
             window = as.integer(chr_covProfile[[x]]$window),
             cumwindow = as.integer(chr_covProfile[[x]]$cumwindow),
             filt_ChIP = as.numeric(filt_ChIP),
             filt_control = as.numeric(filt_control),
             filt_log2ChIPcontrol = as.numeric(filt_log2ChIPcontrol))
}, mc.cores = length(chr_covProfile))

# Combine list of 1 data.frame per chromosome into one data.frame
filt_covProfile <- do.call(rbind, filt_chr_covProfile)

# Write unsmoothed and smoothed chromosome profiles
write.table(covProfile,
            file = paste0(outDir,
                          sampleName, "_", controlName, "_MappedOn_", refbase,
                          "_lowXM_", align, "_sort_norm_binSize", genomeBinName, "_unsmoothed.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(filt_covProfile,
            file = paste0(outDir,
                          sampleName, "_", controlName, "_MappedOn_", refbase,
                          "_lowXM_", align, "_sort_norm_binSize", genomeBinName, "_smoothed.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
