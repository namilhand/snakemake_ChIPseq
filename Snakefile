# Snakemake workflow for mapping ChIP-seq libraries to a reference genome

# Usage (snakemake --cores should reflect available cores):
# conda env create --file environment.yaml --name ChIPseq
# conda activate ChIPseq
# snakemake -p --cores 48
# conda deactivate

# ===== WHAT ELSE AFTER RUNNING THIS PIPELINE? =====
# Use src/bwCompare_TSV.sh to calculate log2 ratio of ChIP-seq signal to control.
# src/genomeBin_bedgraphToTSV.R converts bedgraph form to TSV which can be used for drawing chromosome landscape. TSV separates bedgraph bins that are merged.
# ==================================================

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

# Specify the desired end target file(s)
rule all:
    input:
        expand("qc/fastqc/raw/{sample}",
               sample = sample),
        expand("results/01_trimmed/{sample}_{read}.tr.fastq.gz",
               sample = sample,
               read = [1,2]),
        expand("qc/fastqc/trimmed/{sample}",
               sample = sample),
        expand("results/02_bowtie/{sample}_MappedOn_{refbase}.bam",
                sample = sample,
                refbase = refbase),
        expand("results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.bam",
               sample = sample,
               refbase = refbase),
        expand("results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam",
               sample = sample,
               refbase = refbase),
        expand("results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam.bai",
               sample = sample,
               refbase = refbase),
        expand("results/03_bamCoverage/bw/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_1bp-resolution.bw",
                sample = sample,
                refbase = refbase),
        expand("results/03_bamCoverage/bg/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_1bp-resolution.bedgraph",
               sample = sample,
               refbase = refbase),
        expand("results/03_bamCoverage/bg/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_binSize{genomeBinName}.bedgraph",
               sample = sample,
               refbase = refbase,
               genomeBinName = genomeBinName),
        # expand("results/05_metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}.gz",
        #         sample = ["MNase_suf3-bud", "MNase_WT-bud"],
        #         region = ["genes"],
        #         refpoint = ["TSS", "TES"]),
        # expand("results/05_metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}.tsv",
        #         sample = ["MNase_suf3-bud", "MNase_WT-bud"],
        #         region = ["genes"],
        #         refpoint = ["TSS", "TES"]),

# Run fastqc on paired-end raw data
rule fastqc_raw:
    """Create fastqc report"""
    input:
        read1 = "raw/{sample}_1.fastq.gz",
        read2 = "raw/{sample}_2.fastq.gz"
    output:
        # html = "logs/fastqc/raw/{sample}_fastqc.html",
        # zip  = "logs/fastqc/raw/{sample}_fastqc.zip"
        directory("qc/fastqc/raw/{sample}")
    params:
        " --extract" +
        " --adapters " + str(config["FILTER"]["fastqc"]["adapters"])
    log:
        "logs/fastqc/raw/{sample}.log"
    shell:
        "mkdir -p {output}; fastqc -o {output} {params} {input.read1} {input.read2}"

# Trim off adapters
rule cutadapt:
    """Remove adapters"""
    input:
        read1 = "raw/{sample}_1.fastq.gz",
        read2 = "raw/{sample}_2.fastq.gz"
    output:
        tr1 = "results/01_trimmed/{sample}_1.tr.fastq.gz",
        tr2 = "results/01_trimmed/{sample}_2.tr.fastq.gz",
        qc    = "qc/cutadapt/{sample}_cutadapt.qc.txt"
    params:
        read1_trim=config["FILTER"]["cutadapt"]["read1"], 
        read2_trim=config["FILTER"]["cutadapt"]["read2"], 
        min_overlap=config["FILTER"]["cutadapt"]["minimum-overlap"] 
    log:
        "logs/cutadapt/{sample}_trimmed.log"
    shell:
        r"""
        cutadapt -a {params.read1_trim} -A {params.read2_trim} \
        -O {params.min_overlap} \
        {input.read1} {input.read2} \
        -o {output.tr1} -p {output.tr2} > {output.qc} 2> {log}
        """

# Run fastqc on trimmed data
rule fastqc_trimmed:
    """Create fastqc report"""
    input:
        tr1="results/01_trimmed/{sample}_1.tr.fastq.gz",
        tr2="results/01_trimmed/{sample}_2.tr.fastq.gz"
    output:
        # html = "logs/fastqc/trimmed/{sample}_dedup_trimmed_fastqc.html",
        # zip  = "logs/fastqc/trimmed/{sample}_dedup_trimmed_fastqc.zip"
        directory("qc/fastqc/trimmed/{sample}")
    params:
        " --extract" +
        " --adapters " + str(config["FILTER"]["fastqc"]["adapters"])
    log:
        "logs/fastqc/trimmed/{sample}_trimmed.log"
    shell:
        "mkdir -p {output}; fastqc -o {output} {params} {input.tr1} {input.tr2}"

# Align to reference genome
# Only primary reads with MAPQ > MAPQmaxi are retained
rule bowtie2:
    """Map reads using bowtie2 and filter alignments using samtools"""
    input:
        # fastq = "data/dedup/trimmed/{sample}_dedup_trimmed.fastq.gz",
        tr1="results/01_trimmed/{sample}_1.tr.fastq.gz",
        tr2="results/01_trimmed/{sample}_2.tr.fastq.gz"
    output:
        protected("results/02_bowtie/{sample}_MappedOn_{refbase}.bam")
    params:
        MAPQmaxi = config["MAPPING"]["MAPQmaxi"]
    threads: config["THREADS"]
    log:
        "logs/bowtie2/{sample}_MappedOn_{refbase}.log"
    shell:
        # -F 2308 excludes unmapped reads,
        # as well as secondary and supplementary alignments
        # Exclude alignments with MAPQ < config["MAPPING"]["MAPQmaxi"]
        "(bowtie2 --very-sensitive --no-mixed"
        " --threads {threads}"
        " -x {reference} -1 {input.tr1} -2 {input.tr2}"
        " | samtools view -bh -@ {threads} -F 2308 -q {params.MAPQmaxi} -o {output} - ) 2> {log}"

# Filter alignments for mismatches and extract unique alignments
rule samtools:
    input:"results/02_bowtie/{sample}_MappedOn_{refbase}.bam"
    output:
        protected("results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.bam")
    params:
        sortMemory = config["MAPPING"]["sortMemory"]
    threads: config["THREADS"]
    log:
        "logs/samtools/{sample}_MappedOn_{refbase}_nuclear_sort.log"
    shell:
        r"""
        (samtools view -h -@ {threads} {input} \
        | grep -v -e "ChrM" -e "ChrC" \
        | samtools view -u - \
        | samtools sort -@ {threads} -m {params.sortMemory} -O bam -o {output} -) 2> {log};
        """
        # "| samtools rmdup -s - {output.both}) 2> {log.both}; "
rule markdup:
    output:
        bam="results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam",
        metric="logs/markdup/{sample}_MappedOn_{refbase}_nuclear_sort.md.txt"
        # index="results/02_bowtie2/lowXM/{sample}_MappedOn_{refbase}_lowXM_sort.md.bam.bai"
    input: "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.bam"
    params:
        tmpdir = "tmp"
    shell:
        r"""
        picard MarkDuplicates -I {input} \
        -O {output.bam} \
        -M {output.metric} \
        --TMP_DIR {params.tmpdir} \
        --REMOVE_DUPLICATES true;
        """
rule postmapping:
    """bam.bai samtools flagstat idxstats"""
    input:
        "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam"
    output:
        "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam.bai"
    log:
        flagstat   = "qc/samtools/stats/{sample}_MappedOn_{refbase}_nuclear_sort_md_flagstat.log",
        idxstats   = "qc/samtools/stats/{sample}_MappedOn_{refbase}_nuclear_sort_md_idxstats.log"
    shell:
        """
        samtools index    {input}
        samtools flagstat {input} > {log.flagstat}
        samtools idxstats {input} > {log.idxstats}
        """
rule calc_coverage:
    """Calculate library-size-normalized coverage"""
    input:
        BAM   = "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam",
        BAMidx   = "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam.bai"
    output:
        BW   = "results/03_bamCoverage/bw/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_1bp-resolution.bw",
        BG   = "results/03_bamCoverage/bg/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_1bp-resolution.bedgraph"
    params:
        # extendReads option for paired end sequencing:
        # properly paired -> extend to match the fragment size
        # not paired or not properly paired -> treated as singled-end reads. extend to specified length.
        normalizeUsing         = config["COVERAGE"]["normalizeUsing"],
        ignoreForNormalization = config["COVERAGE"]["ignoreForNormalization"],
        extendReads            = config["COVERAGE"]["extendReads"],
        binSize                = config["COVERAGE"]["binSize"]
    log:
        "logs/bamCoverage/{sample}_MappedOn_{refbase}_nuclear_both_sort_md_norm.log"
    threads: config["THREADS"]  
    shell:
        "(bamCoverage -b {input.BAM} -o {output.BW}"
        " --normalizeUsing {params.normalizeUsing}"
        " --ignoreForNormalization {params.ignoreForNormalization}"
        " --exactScaling"
        " --extendReads {params.extendReads}"
        " --binSize {params.binSize} -p {threads}; "
        "bamCoverage -b {input.BAM} -o {output.BG} -of bedgraph"
        " --normalizeUsing {params.normalizeUsing}"
        " --ignoreForNormalization {params.ignoreForNormalization}"
        " --exactScaling"
        " --extendReads {params.extendReads}"
        " --binSize {params.binSize} -p {threads}) 2> {log}"
rule calc_coverage_genome:
    """Calculate library-size-normalized coverage in adjacent windows"""
    input:
        BAM   = "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam",
        BAMidx   = "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam.bai"
    output:
        "results/03_bamCoverage/bg/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_binSize{genomeBinName}.bedgraph"
    params:
        normalizeUsing         = config["COVERAGE"]["normalizeUsing"],
        ignoreForNormalization = config["COVERAGE"]["ignoreForNormalization"],
        extendReads            = config["COVERAGE"]["extendReads"],
        genomeBinSize          = config["COVERAGE"]["genomeBinSize"]
    log:
        "logs/bamCoverage/{sample}_MappedOn_{refbase}_nuclear_both_sort_norm_binSize{genomeBinName}.log"
    threads: config["THREADS"]  
    shell:
        "(bamCoverage -b {input.BAM} -o {output} -of bedgraph"
        " --normalizeUsing {params.normalizeUsing}"
        " --ignoreForNormalization {params.ignoreForNormalization}"
        " --exactScaling"
        " --extendReads {params.extendReads}"
        " --binSize {params.genomeBinSize} -p {threads}) 2> {log}"
rule metaprofile:
    output:
        gzip = "results/05_metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}.gz",
        tab = "results/05_metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}.tab",
        bed = "results/05_metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}.bed"
    input:
        "results/04_log2ChIP/{sample}_log2ChIP_binSize1bp.bw"
    params:
        refpoint="{refpoint}",
        region = config["METAPROFILE"]["region"],
        binSize = config["METAPROFILE"]["binSize"]
    threads: config["THREADS"] 
    log:
        "logs/metaprofile/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}.log"
    shell:
        r"""
        computeMatrix reference-point \
            --referencePoint {params.refpoint} \
            -b 2000 -a 2000 \
            -R {params.region} \
            -S {input} \
            --binSize {params.binSize} \
            --sortRegions "keep" \
            --sortUsing "mean" \
            --samplesLabel {wildcards.sample} \
            --nanAfterEnd \
            -p {threads} \
            --outFileName {output.gzip} \
            --outFileNameMatrix {output.tab} \
            --outFileSortedRegions {output.bed} &> {log}
        """
rule avgBins:
    output: "results/05_metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}.tsv"
    input: "results/05_metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}.tab"
    params:
        upstream = config["METAPROFILE"]["upstream"],
        downstream = config["METAPROFILE"]["downstream"],
        binSize = config["METAPROFILE"]["binSize"]
    shell:
        r"""
        Rscript src/averageComputeMatrixTab.R {input} \
        {params.upstream} {params.downstream} \
        {params.binSize} \
        {output}
        """