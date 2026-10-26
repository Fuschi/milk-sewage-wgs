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

# Detect chimerism in DAS Tool bins.
rule gunc_bins:
    input:
        bins="data/refinement/dastool/{mode}/{unit}/{unit}_DASTool_bins",
        assignments="data/refinement/dastool/{mode}/{unit}/{unit}_DASTool_contig2bin.tsv",
        db=config["references"]["gunc"]
    output: "data/bin_quality/gunc/{mode}/{unit}/GUNC.progenomes_3.maxCSS_level.tsv"
    params: outdir="data/bin_quality/gunc/{mode}/{unit}"
    wildcard_constraints: mode="sample|coassembly"
    threads: 16
    resources: mem_mb=64000, runtime=480
    container: config["containers"]["gunc"]
    log: "logs/bin_quality/gunc/{mode}/{unit}.log"
    benchmark: "benchmarks/bin_quality/gunc/{mode}/{unit}.tsv"
    shell:
        """
        if [ -s {input.assignments} ]; then
            gunc run --input_dir {input.bins} --file_suffix .fa --db_file {input.db} --out_dir {params.outdir} --threads {threads} --detailed_output > {log} 2>&1
        else
            echo "No DAS Tool bins to assess." > {log}
            : > {output}
        fi
        """
