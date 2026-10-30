# Profile dereplicated MAGs with CoverM and ABRicate.

rule coverm_mags:
    input:
        drep="data/dereplication/drep",
        r1="data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz",
        r2="data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz"
    output:
        "data/mags_profiling/coverm/{sample}.tsv"
    params:
        genomes="data/dereplication/drep/dereplicated_genomes",
        methods=(
            "relative_abundance mean trimmed_mean covered_fraction "
            "covered_bases variance rpkm tpm"
        )
    threads: 16
    resources:
        mem_mb=32000,
        runtime=720
    container:
        config["containers"]["coverm"]
    log:
        "logs/mags_profiling/coverm/{sample}.log"
    benchmark:
        "benchmarks/mags_profiling/coverm/{sample}.tsv"
    shell:
        """
        coverm genome \
            --genome-fasta-directory {params.genomes} \
            --genome-fasta-extension fa \
            --coupled {input.r1} {input.r2} \
            --methods {params.methods} \
            --threads {threads} \
            --output-file {output} \
            > {log} 2>&1
        """

rule abricate_mags:
    input:
        drep="data/dereplication/drep"
    output:
        "tables/mags_profiling/abricate_{database}.tsv"
    params:
        genomes="data/dereplication/drep/dereplicated_genomes",
        min_identity=90,
        min_coverage=60
    wildcard_constraints:
        database="resfinder_2.6.0|vfdb_2026-10-26"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=120
    container:
        config["containers"]["abricate"]
    log:
        "logs/mags_profiling/abricate/{database}.log"
    benchmark:
        "benchmarks/mags_profiling/abricate/{database}.tsv"
    shell:
        """
        abricate \
            --db {wildcards.database} \
            --minid {params.min_identity} \
            --mincov {params.min_coverage} \
            --threads {threads} \
            --nopath \
            {params.genomes}/*.fa \
            > {output} 2> {log}
        """

# Combine CoverM results across samples.

rule coverm_mags_summary:
    input:
        coverm=expand(
            "data/mags_profiling/coverm/{sample}.tsv",
            sample=SAMPLES
        )
    output:
        summary="tables/mags_profiling/mags_abundance.tsv"
    resources:
        mem_mb=4000,
        runtime=30
    container:
        config["containers"]["r"]
    log:
        "logs/mags_profiling/coverm_summary.log"
    benchmark:
        "benchmarks/mags_profiling/coverm_summary.tsv"
    script:
        "../scripts/aggregate_coverm_mags.R"
