# Snakemake workflow for mapping ChIP-seq libraries to a reference genome

# Chromosome sizes file below ("data/index/genome.fa.sizes") must exist
# before running snakefile
# e.g., in "data/index/" run:
# samtools faidx genome.fa; cut -f1,2 genome.fa.fai > genome.fa.sizes

# Usage (snakemake --cores should reflect available cores):
# conda env create --file environment.yaml --name ChIPseq_mapping
# conda activate ChIPseq_mapping
# snakemake -p --cores 48
# conda deactivate

import pandas as pd
import os

# To make the per_base_coverage rule work with a shell script invoked using the "shell" directive,
# we need to determine the base path of Snakefile since we expect the scripts directory to be there as well
SRCDIR = srcdir("")

# Specify config file parameters
configfile: "config.yaml"
sample        = config["SAMPLES"]
reference     = config["MAPPING"]["reference"]
refbase       = os.path.basename(reference)
genomeBinName = config["COVERAGE"]["genomeBinName"]

# Determine bam index format (bai or csi) based on chromosome sizes
# Genomes with chromosomes longer than ~500 Mb (e.g., in wheat) require a csi index
# E.g., in axolotl: https://sourceforge.net/p/samtools/mailman/message/36249039/
chrSizes = pd.read_csv("data/index/" + refbase + ".fa.sizes",
                       header = None, sep = "\t")
smallChrs = 0
for x in chrSizes[1]:
    if x < 5e+08:
        smallChrs = smallChrs + 1

if smallChrs < len(chrSizes[1]):
    bamidx = "csi"
else:
    bamidx = "bai"

# Specify the desired end target file(s)
rule all:
    input:
        expand("logs/fastqc/raw/{sample}_fastqc.html",
               sample = sample),
        expand("data/dedup/{sample}_dedup.fastq.gz",
               sample = sample),
        expand("data/dedup/trimmed/{sample}_dedup_trimmed.fastq.gz",
               sample = sample),
        expand("logs/fastqc/trimmed/{sample}_dedup_trimmed_fastqc.html",
               sample = sample),
        expand("mapped/{sample}_MappedOn_{refbase}.bam",
               sample = sample,
               refbase = refbase),
        expand("mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam",
               sample = sample,
               refbase = refbase),
        expand("mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam",
               sample = sample,
               refbase = refbase),
        expand("mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam.{bamidx}",
               sample = sample,
               refbase = refbase,
               bamidx = bamidx),
        expand("mapped/unique/bw/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.bw",
               sample = sample,
               refbase = refbase),
        expand("mapped/unique/bg/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.bedgraph",
               sample = sample,
               refbase = refbase),
        expand("mapped/unique/bg/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.bedgraph",
               sample = sample,
               refbase = refbase,
               genomeBinName = genomeBinName),
        expand("mapped/unique/tsv/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.tsv",
               sample = sample,
               refbase = refbase,
               genomeBinName = genomeBinName),
        expand("mapped/unique/pb/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.perbase",
               sample = sample,
               refbase = refbase),
        expand("mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam.{bamidx}",
               sample = sample,
               refbase = refbase,
               bamidx = bamidx),
        expand("mapped/both/bw/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.bw",
               sample = sample,
               refbase = refbase),
        expand("mapped/both/bg/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.bedgraph",
               sample = sample,
               refbase = refbase),
        expand("mapped/both/bg/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.bedgraph",
               sample = sample,
               refbase = refbase,
               genomeBinName = genomeBinName),
        expand("mapped/both/tsv/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.tsv",
               sample = sample,
               refbase = refbase,
               genomeBinName = genomeBinName),
        expand("mapped/both/pb/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.perbase",
               sample = sample,
               refbase = refbase)

# Run fastqc on single-end raw data
rule fastqc_single_raw:
    """Create fastqc report"""
    input:
        "data/{sample}.fastq.gz"
    output:
        html = "logs/fastqc/raw/{sample}_fastqc.html",
        zip  = "logs/fastqc/raw/{sample}_fastqc.zip"
    params:
        " --extract" +
        " --adapters " + str(config["FILTER"]["fastqc"]["adapters"])
    log:
        "logs/fastqc/raw/{sample}.log"
    wrapper:
        "0.27.1/bio/fastqc"

# Deduplicate single-end reads
rule dedupe_single:
    """Remove duplicate single-end reads"""
    input:
        "data/{sample}.fastq.gz"
    output:
        "data/dedup/{sample}_dedup.fastq.gz"
    threads: config["THREADS"]
    params:
        memory = config["MEMORY"]
    log:
        "logs/dedup/{sample}_dedup.log"
    shell:
        "(dedupe.sh -Xmx{params.memory} in={input} out={output}"
        " threads={threads} ac=f) 2> {log}"

# Trim off adapters
rule cutadapt:
    """Remove adapters"""
    input:
        "data/dedup/{sample}_dedup.fastq.gz"
    output:
        fastq = "data/dedup/trimmed/{sample}_dedup_trimmed.fastq.gz",
        qc    = "data/dedup/trimmed/{sample}_dedup_trimmed.qc.txt"
    params:
        " -u " + str(config["FILTER"]["cutadapt"]["5prime_cut"]) +
        " -u " + str(config["FILTER"]["cutadapt"]["3prime_cut"]) +
        " -a " +     config["FILTER"]["cutadapt"]["adapter"] +
        " -O " + str(config["FILTER"]["cutadapt"]["minimum-overlap"]) +
        " -q " + str(config["FILTER"]["cutadapt"]["quality-filter"]) +
        " -m " + str(config["FILTER"]["cutadapt"]["minimum-length"]) +
        " -M " + str(config["FILTER"]["cutadapt"]["maximum-length"]) +
        " --cores=0"
    log:
        "logs/cutadapt/{sample}_dedup_trimmed.log"
    wrapper:
        "0.27.1/bio/cutadapt/se"

# Run fastqc on single-end trimmed data
rule fastqc_single_trimmed:
    """Create fastqc report"""
    input:
        "data/dedup/trimmed/{sample}_dedup_trimmed.fastq.gz"
    output:
        html = "logs/fastqc/trimmed/{sample}_dedup_trimmed_fastqc.html",
        zip  = "logs/fastqc/trimmed/{sample}_dedup_trimmed_fastqc.zip"
    params:
        " --extract" +
        " --adapters " + str(config["FILTER"]["fastqc"]["adapters"])
    log:
        "logs/fastqc/trimmed/{sample}_dedup_trimmed.log"
    wrapper:
        "0.27.1/bio/fastqc"

# Align to reference genome
rule bowtie2:
    """Map reads using bowtie2 and filter alignments using samtools"""
    input:
        fastq = "data/dedup/trimmed/{sample}_dedup_trimmed.fastq.gz",
    output:
        protected("mapped/{sample}_MappedOn_{refbase}.bam")
    params:
        alignments = config["MAPPING"]["alignments"],
        MAPQmaxi = config["MAPPING"]["MAPQmaxi"]
    threads: config["THREADS"]
    log:
        "logs/bowtie2/{sample}_MappedOn_{refbase}_sort.log"
    shell:
        # -F 2308 excludes unmapped reads,
        # as well as secondary and supplementary alignments
        # Exclude alignments with MAPQ < config["MAPPING"]["MAPQmaxi"]
        "(bowtie2 --very-sensitive"
        " --threads {threads} -k {params.alignments}"
        " -x {reference} -U {input.fastq} "
        "| samtools view -bh -@ {threads} -F 2308 -q {params.MAPQmaxi} -o {output} - ) 2> {log}"

# Filter alignments for mismatches and extract unique alignments
rule samtools:
    input:
        "mapped/{sample}_MappedOn_{refbase}.bam"
    output:
        both   = protected("mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam"),
        unique = protected("mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam")
    params:
        sortMemory = config["MAPPING"]["sortMemory"],
        MAPQunique = config["MAPPING"]["MAPQunique"]
    threads: config["THREADS"]
    log:
        both   = "logs/samtools/{sample}_MappedOn_{refbase}_lowXM_both_sort.log",
        unique = "logs/samtools/{sample}_MappedOn_{refbase}_lowXM_unique_sort.log"
    shell:
        # Allow a maximum of 2 mismatches
        # ([^0-9] matches characters not in the range of 0 to 9)
        # http://seqanswers.com/forums/showthread.php?t=19729
        "(samtools view -h {input} " 
        "| grep -e '^@' -e 'XM:i:[012][^0-9]' "
        "| samtools view -u - "
        "| samtools sort -@ {threads} -m {params.sortMemory} - "
        "| samtools rmdup -s - {output.both}) 2> {log.both}; "
        # Extract unique alignments, excluding alignments with MAPQ scores < config["MAPPING"]["MAPQunique"]
        # http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
        # https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/
        "(samtools view -h -q {params.MAPQunique} {input} "
        "| grep -e '^@' -e 'XM:i:[012][^0-9]' "
        "| samtools view -u - "
        "| samtools sort -@ {threads} -m {params.sortMemory} - "
        "| samtools rmdup -s - {output.unique}) 2> {log.unique}"

# Postmapping steps:
# Index BAM files (index format [bai or csi] depends on chromosome sizes)
# Generate samtools flagstat and idxstats
# Calculate library-size-normalized coverage
if bamidx == "bai":
    rule postmapping:
        """bam.bai samtools flagstat idxstats"""
        input:
            uniqueBAM = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam",
            bothBAM   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam"
        output:
            uniqueBAM = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam.bai",
            bothBAM   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam.bai"
        log:
            uniqueflagstat = "logs/samtools/stats/{sample}_MappedOn_{refbase}_lowXM_unique_sort_flagstat.log",
            bothflagstat   = "logs/samtools/stats/{sample}_MappedOn_{refbase}_lowXM_both_sort_flagstat.log",
            uniqueidxstats = "logs/samtools/stats/{sample}_MappedOn_{refbase}_lowXM_unique_sort_idxstats.log",
            bothidxstats   = "logs/samtools/stats/{sample}_MappedOn_{refbase}_lowXM_both_sort_idxstats.log"
        shell:
            """
            samtools index    {input.uniqueBAM}
            samtools flagstat {input.uniqueBAM} > {log.uniqueflagstat}
            samtools idxstats {input.uniqueBAM} > {log.uniqueidxstats}
            samtools index    {input.bothBAM}
            samtools flagstat {input.bothBAM} > {log.bothflagstat}
            samtools idxstats {input.bothBAM} > {log.bothidxstats}
            """
    rule calc_coverage:
        """Calculate library-size-normalized coverage"""
        input:
            uniqueBAM = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam",
            bothBAM   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam",
            uniqueBAMidx = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam.bai",
            bothBAMidx   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam.bai"
        output:
            uniqueBW = "mapped/unique/bw/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.bw",
            bothBW   = "mapped/both/bw/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.bw",
            uniqueBG = "mapped/unique/bg/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.bedgraph",
            bothBG   = "mapped/both/bg/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.bedgraph"
        params:
            normalizeUsing         = config["COVERAGE"]["normalizeUsing"],
            ignoreForNormalization = config["COVERAGE"]["ignoreForNormalization"],
            extendReads            = config["COVERAGE"]["extendReads"],
            binSize                = config["COVERAGE"]["binSize"]
        log:
            unique = "logs/bamCoverage/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.log",
            both   = "logs/bamCoverage/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.log"
        threads: config["THREADS"]  
        shell:
            "(bamCoverage -b {input.uniqueBAM} -o {output.uniqueBW}"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.binSize} -p {threads}; "
            "bamCoverage -b {input.uniqueBAM} -o {output.uniqueBG} -of bedgraph"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.binSize} -p {threads}) 2> {log.unique}; "
            "(bamCoverage -b {input.bothBAM} -o {output.bothBW}"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.binSize} -p {threads}; "
            "bamCoverage -b {input.bothBAM} -o {output.bothBG} -of bedgraph"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.binSize} -p {threads}) 2> {log.both}"
    rule calc_coverage_genome:
        """Calculate library-size-normalized coverage in adjacent windows"""
        input:
            uniqueBAM = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam",
            bothBAM   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam",
            uniqueBAMidx = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam.bai",
            bothBAMidx   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam.bai"
        output:
            uniqueBGgenome = "mapped/unique/bg/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.bedgraph",
            bothBGgenome   = "mapped/both/bg/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.bedgraph"
        params:
            normalizeUsing         = config["COVERAGE"]["normalizeUsing"],
            ignoreForNormalization = config["COVERAGE"]["ignoreForNormalization"],
            extendReads            = config["COVERAGE"]["extendReads"],
            genomeBinSize          = config["COVERAGE"]["genomeBinSize"]
        log:
            unique = "logs/bamCoverage/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.log",
            both   = "logs/bamCoverage/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.log"
        threads: config["THREADS"]  
        shell:
            "(bamCoverage -b {input.uniqueBAM} -o {output.uniqueBGgenome} -of bedgraph"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.genomeBinSize} -p {threads}) 2> {log.unique}; "
            "(bamCoverage -b {input.bothBAM} -o {output.bothBGgenome} -of bedgraph"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.genomeBinSize} -p {threads}) 2> {log.both}"
else:
    rule postmapping:
        """bam.csi samtools flagstat idxstats"""
        input:
            uniqueBAM = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam",
            bothBAM   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam"
        output:
            uniqueBAM = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam.csi",
            bothBAM   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam.csi"
        log:
            uniqueflagstat = "logs/samtools/stats/{sample}_MappedOn_{refbase}_lowXM_unique_sort_flagstat.log",
            bothflagstat   = "logs/samtools/stats/{sample}_MappedOn_{refbase}_lowXM_both_sort_flagstat.log",
            uniqueidxstats = "logs/samtools/stats/{sample}_MappedOn_{refbase}_lowXM_unique_sort_idxstats.log",
            bothidxstats   = "logs/samtools/stats/{sample}_MappedOn_{refbase}_lowXM_both_sort_idxstats.log"
        shell:
            """
            samtools index -c -m 14 {input.uniqueBAM}
            samtools flagstat       {input.uniqueBAM} > {log.uniqueflagstat}
            samtools idxstats       {input.uniqueBAM} > {log.uniqueidxstats}
            samtools index -c -m 14 {input.bothBAM}
            samtools flagstat       {input.bothBAM} > {log.bothflagstat}
            samtools idxstats       {input.bothBAM} > {log.bothidxstats}
            """
    rule calc_coverage:
        """Calculate library-size-normalized coverage"""
        input:
            uniqueBAM = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam",
            bothBAM   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam",
            uniqueBAMidx = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam.csi",
            bothBAMidx   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam.csi"
        output:
            uniqueBW = "mapped/unique/bw/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.bw",
            bothBW   = "mapped/both/bw/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.bw",
            uniqueBG = "mapped/unique/bg/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.bedgraph",
            bothBG   = "mapped/both/bg/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.bedgraph"
        params:
            normalizeUsing         = config["COVERAGE"]["normalizeUsing"],
            ignoreForNormalization = config["COVERAGE"]["ignoreForNormalization"],
            extendReads            = config["COVERAGE"]["extendReads"],
            binSize                = config["COVERAGE"]["binSize"]
        log:
            unique = "logs/bamCoverage/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.log",
            both   = "logs/bamCoverage/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.log"
        threads: config["THREADS"]  
        shell:
            "(bamCoverage -b {input.uniqueBAM} -o {output.uniqueBW}"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.binSize} -p {threads}; "
            "bamCoverage -b {input.uniqueBAM} -o {output.uniqueBG} -of bedgraph"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.binSize} -p {threads}) 2> {log.unique}; "
            "(bamCoverage -b {input.bothBAM} -o {output.bothBW}"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.binSize} -p {threads}; "
            "bamCoverage -b {input.bothBAM} -o {output.bothBG} -of bedgraph"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.binSize} -p {threads}) 2> {log.both}"
    rule calc_coverage_genome:
        """Calculate library-size-normalized coverage in adjacent windows"""
        input:
            uniqueBAM = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam",
            bothBAM   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam",
            uniqueBAMidx = "mapped/unique/{sample}_MappedOn_{refbase}_lowXM_unique_sort.bam.csi",
            bothBAMidx   = "mapped/both/{sample}_MappedOn_{refbase}_lowXM_both_sort.bam.csi"
        output:
            uniqueBGgenome = "mapped/unique/bg/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.bedgraph",
            bothBGgenome   = "mapped/both/bg/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.bedgraph"
        params:
            normalizeUsing         = config["COVERAGE"]["normalizeUsing"],
            ignoreForNormalization = config["COVERAGE"]["ignoreForNormalization"],
            extendReads            = config["COVERAGE"]["extendReads"],
            genomeBinSize          = config["COVERAGE"]["genomeBinSize"]
        log:
            unique = "logs/bamCoverage/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.log",
            both   = "logs/bamCoverage/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.log"
        threads: config["THREADS"]  
        shell:
            "(bamCoverage -b {input.uniqueBAM} -o {output.uniqueBGgenome} -of bedgraph"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.genomeBinSize} -p {threads}) 2> {log.unique}; "
            "(bamCoverage -b {input.bothBAM} -o {output.bothBGgenome} -of bedgraph"
            " --normalizeUsing {params.normalizeUsing}"
            " --ignoreForNormalization {params.ignoreForNormalization}"
            " --exactScaling"
            " --extendReads {params.extendReads}"
            " --binSize {params.genomeBinSize} -p {threads}) 2> {log.both}"

# Use R script genomeBin_bedgraphToTSV.R to convert *{genomeBinName}.bedgraph files into TSV files
# These TSV files can be imported into R for calculating and plotting log2(ChIP/control) chromosome-scale profiles
rule bedgraphToTSV:
    """Convert *{genomeBinName}.bedgraph files into TSV files"""
    input:
        uniqueBGgenome = "mapped/unique/bg/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.bedgraph",
        bothBGgenome   = "mapped/both/bg/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.bedgraph"
    output:
        uniqueTSVgenome = "mapped/unique/tsv/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.tsv",
        bothTSVgenome   = "mapped/both/tsv/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.tsv"
    params:
        genomeBinSize = config["COVERAGE"]["genomeBinSize"]
    log:
        unique = "logs/genomeBin_bedgraphToTSV/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.log",
        both   = "logs/genomeBin_bedgraphToTSV/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.log"
    threads: config["THREADS"]
    shell:
        "(scripts/genomeBin_bedgraphToTSV.R"
        " {wildcards.sample}"
        " {refbase}"
        " unique"
        " {params.genomeBinSize}) 2> {log.unique}; "
        "(scripts/genomeBin_bedgraphToTSV.R"
        " {wildcards.sample}"
        " {refbase}"
        " both"
        " {params.genomeBinSize}) 2> {log.both}"

# Convert bedgraph to per-base 1-based coverage file
rule per_base_coverage:
    """Convert bedgraph to per-base 1-based coverage file"""
    input:
        uniqueBG = "mapped/unique/bg/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.bedgraph",
        bothBG   = "mapped/both/bg/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.bedgraph"
    output:
        uniquePB = "mapped/unique/pb/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.perbase",
        bothPB   = "mapped/both/pb/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.perbase"
    log:
        unique = "logs/perBaseCoverage/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_pb.log",
        both = "logs/perBaseCoverage/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_pb.log"
    shell:
        "(bash scripts/perbase_1based_coverage.sh {input.bothBG} {output.bothPB} ) 2> {log.both}; "
        "(bash scripts/perbase_1based_coverage.sh {input.uniqueBG} {output.uniquePB} ) 2> {log.unique}"
