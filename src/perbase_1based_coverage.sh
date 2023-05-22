#!/bin/bash

# Convert bedgraph to per-base 1-based coverage file

# Usage from within Snakefile
#rule per_base_coverage:
#    """Convert bedgraph to per-base 1-based coverage file"""
#    input:
#        "mapped/bg/{sample}_MappedOn_{refbase}_{mode}_norm.bedgraph"
#    output:
#        "mapped/pb/{sample}_MappedOn_{refbase}_{mode}_norm.perbase"
#    log:
#        "logs/perBaseCoverage/{sample}_MappedOn_{refbase}_{mode}_norm_pb.log"
#    shell:
#        "bash {SRCDIR}/scripts/perbase_1based_coverage.sh {input} {output}"

perl -alne 'for ($F[1]+1..$F[2]) { print "$F[0]\t$_\t$F[3]" }' $1 > $2
