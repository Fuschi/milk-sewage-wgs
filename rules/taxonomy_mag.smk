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
