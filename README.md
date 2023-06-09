# Snakemake pipeline for general ChIP-seq analysis

The snakefile includes ChIP-seq analysis processes from fastqc to bamCoverage, resulting in genomewide coverage of ChIP-seq data (bigwig).

Normalization of ChIP-seq reads by control data (such as ChIP-seq input) is not included in this pipeline as the process is tricky to be generallized. Instead, you can use src/bwCompare_TSV.sh for normalizing ChIP-seq data. 

# TL;DR
1. conda env create --name ChIPseq --file environment.yaml
2. conda activate ChIPseq
3. Edit config.yaml
4. snakemake -p --cores <threads>
5. Apply the step 2-4 to both control and ChIP-seq data.
6. Edit and run src/bwCompare_TSV.sh
7. If you want to draw a metaprofile of a set of region (deeptools:computeMatrix), uncomment rule "metaprofile" and "avgBins". Edit Snakefile input and configure.yaml accordingly, and run the Snakemake again.
