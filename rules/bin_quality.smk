# Estimate completeness and contamination of DAS Tool bins.
rule checkm2_bins:
    input:
        bins="data/refinement/dastool/{mode}/{unit}/{unit}_DASTool_bins",
        assignments="data/refinement/dastool/{mode}/{unit}/{unit}_DASTool_contig2bin.tsv",
        db=config["references"]["checkm2"]
    output: "data/bin_quality/checkm2/{mode}/{unit}/quality_report.tsv"
    params: outdir="data/bin_quality/checkm2/{mode}/{unit}"
    wildcard_constraints: mode="sample|coassembly"
    threads: 16
    resources: mem_mb=64000, runtime=480
    container: config["containers"]["checkm2"]
    log: "logs/bin_quality/checkm2/{mode}/{unit}.log"
    benchmark: "benchmarks/bin_quality/checkm2/{mode}/{unit}.tsv"
    shell:
        """
        if [ -s {input.assignments} ]; then
            checkm2 predict --input {input.bins} --extension fa --output-directory {params.outdir} --database_path {input.db} --threads {threads} > {log} 2>&1
        else
            echo "No DAS Tool bins to assess." > {log}
            : > {output}
        fi
        """

# Prepare the bin quality summary and dRep inputs.
rule prepare_drep_inputs:
    input:
        bins_sample=expand("data/refinement/dastool/sample/{sample}/{sample}_DASTool_bins", sample=SAMPLES),
        bins_coassembly=expand("data/refinement/dastool/coassembly/{biome}/{biome}_DASTool_bins", biome=BIOMES),
        checkm2_sample=expand("data/bin_quality/checkm2/sample/{sample}/quality_report.tsv", sample=SAMPLES),
        checkm2_coassembly=expand("data/bin_quality/checkm2/coassembly/{biome}/quality_report.tsv", biome=BIOMES),
        gunc_sample=expand("data/bin_quality/gunc/sample/{sample}/GUNC.progenomes_3.maxCSS_level.tsv", sample=SAMPLES),
        gunc_coassembly=expand("data/bin_quality/gunc/coassembly/{biome}/GUNC.progenomes_3.maxCSS_level.tsv", biome=BIOMES)
    output:
        summary="tables/bin_quality_summary.tsv",
        genomes="data/dereplication/inputs/genomes.txt",
        genome_info="data/dereplication/inputs/genome_info.csv"
    container: config["containers"]["r"]
    script:
        "../scripts/prepare_drep_inputs.R"
