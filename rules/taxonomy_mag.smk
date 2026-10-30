# Classify dereplicated MAGs with GTDB-Tk.
rule gtdbtk:
    input:
        drep="data/drep",
        db=config["references"]["gtdbtk"]
    output: directory("data/gtdbtk")
    params: genomes="data/drep/dereplicated_genomes"
    threads: 64
    resources:
        mem_mb=1000000,
        runtime=1440,
        slurm_partition="m8",
        slurm_extra="--nodelist=mtx30 --reservation=prj-trigger"
    container: config["containers"]["gtdbtk"]
    log: "logs/taxonomy/gtdbtk.log"
    benchmark: "benchmarks/taxonomy/gtdbtk.tsv"
    shell:
        """
        gtdbtk classify_wf \
            --genome_dir {params.genomes} \
            --out_dir {output} \
            --extension fa \
            --cpus {threads} \
            --full_tree \
            > {log} 2>&1
        """

# Combine bacterial and archaeal GTDB-Tk classifications.
rule gtdbtk_summary:
    input:
        ar53="data/gtdbtk/classify/gtdbtk.ar53.summary.tsv",
        bac120="data/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
    output: summary="results/gtdbtk_summary.tsv"
    resources:
        mem_mb=4000,
        runtime=30
    container: config["containers"]["r"]
    benchmark: "benchmarks/taxonomy/gtdbtk_summary.tsv"
    script: "../scripts/aggregate_gtdbtk.R"

# Infer NCBI taxonomy from GTDB-Tk classifications.
rule gtdbtk_ncbi_taxonomy:
    input:
        gtdbtk="data/gtdbtk",
        ar53=config["references"]["gtdb_ar53_metadata"],
        bac120=config["references"]["gtdb_bac120_metadata"]
    output:
        temp("data/gtdbtk_ncbi_raw.tsv")
    threads: 1
    resources:
        mem_mb=64000,
        runtime=240
    container:
        config["containers"]["gtdbtk"]
    log:
        "logs/taxonomy/gtdbtk_ncbi_taxonomy.log"
    benchmark:
        "benchmarks/taxonomy/gtdbtk_ncbi_taxonomy.tsv"
    shell:
        """
        gtdb_to_ncbi_majority_vote.py \
            --gtdbtk_output_dir {input.gtdbtk:q} \
            --ar53_metadata_file {input.ar53:q} \
            --bac120_metadata_file {input.bac120:q} \
            --output_file {output:q} \
            > {log:q} 2>&1
        """

# Prepare the NCBI taxonomy table for downstream analysis.
rule gtdbtk_ncbi_summary:
    input:
        "data/gtdbtk_ncbi_raw.tsv"
    output:
        "results/ncbi_taxonomy.tsv"
    resources:
        mem_mb=2000,
        runtime=10
    container:
        config["containers"]["r"]
    benchmark:
        "benchmarks/taxonomy/gtdbtk_ncbi_summary.tsv"
    script:
        "../scripts/aggregate_ncbi_taxonomy.R"

