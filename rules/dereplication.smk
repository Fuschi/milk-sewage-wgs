# Dereplicate selected bins and choose representative genomes.
rule drep:
    input:
        genomes="data/drep_inputs/genomes.txt",
        genome_info="data/drep_inputs/genome_info.csv"
    output: directory("data/drep")
    threads: 32
    resources:
        mem_mb=128000,
        runtime=1440
    container: config["containers"]["drep"]
    log: "logs/dereplication/drep.log"
    benchmark: "benchmarks/dereplication/drep.tsv"
    shell:
        """
        dRep dereplicate {output} \
            -g {input.genomes} \
            --genomeInfo {input.genome_info} \
            -p {threads} \
            -comp 50 \
            -con 10 \
            -sa 0.95 \
            > {log} 2>&1
        """
