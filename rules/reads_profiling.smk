rule profile_amr_reads:
    input:
        r1="data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz",
        r2="data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz",
        s1="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz",
        s2="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz",
        db=expand(config["references"]["resfinder_kma"] + ".{suffix}", suffix=["comp.b", "length.b", "name", "seq.b"])
    output:
        res="data/reads_profiling/amr/resfinder/{sample}/{sample}.res",
        mapstat="data/reads_profiling/amr/resfinder/{sample}/{sample}.mapstat"
    params:
        db=config["references"]["resfinder_kma"],
        prefix="data/reads_profiling/amr/resfinder/{sample}/{sample}"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=120
    container: config["containers"]["kma"]
    log: "logs/reads_profiling/amr/resfinder/{sample}.log"
    benchmark: "benchmarks/reads_profiling/amr/resfinder/{sample}.tsv"
    shell:
        """
        kma -ipe {input.r1} {input.r2} -i {input.s1} {input.s2} \
            -t_db {params.db} -o {params.prefix} -e 1.0 -apm p -cge -1t1 -t {threads} -ef > {log} 2>&1
        """


rule profile_virulence_reads:
    input:
        r1="data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz",
        r2="data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz",
        s1="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz",
        s2="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz",
        db=expand(config["references"]["vfdb_kma"] + ".{suffix}", suffix=["comp.b", "length.b", "name", "seq.b"])
    output:
        res="data/reads_profiling/virulence/vfdb/{sample}/{sample}.res",
        mapstat="data/reads_profiling/virulence/vfdb/{sample}/{sample}.mapstat"
    params:
        db=config["references"]["vfdb_kma"],
        prefix="data/reads_profiling/virulence/vfdb/{sample}/{sample}"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=120
    container: config["containers"]["kma"]
    log: "logs/reads_profiling/virulence/vfdb/{sample}.log"
    benchmark: "benchmarks/reads_profiling/virulence/vfdb/{sample}.tsv"
    shell:
        """
        kma -ipe {input.r1} {input.r2} -i {input.s1} {input.s2} \
            -t_db {params.db} -o {params.prefix} -e 1.0 -apm p -cge -1t1 -t {threads} -ef > {log} 2>&1
        """


rule quantify_amr_reads:
    input:
        res="data/reads_profiling/amr/resfinder/{sample}/{sample}.res",
        mapstat="data/reads_profiling/amr/resfinder/{sample}/{sample}.mapstat",
        script="scripts/quantify_kma_reads.sh"
    output:
        all_refs="data/reads_profiling/amr/resfinder/{sample}/{sample}.quantification.all.tsv",
        filtered="data/reads_profiling/amr/resfinder/{sample}/{sample}.quantification.filtered.tsv"
    params:
        prefix="data/reads_profiling/amr/resfinder/{sample}/{sample}.quantification",
        min_identity=90,
        min_coverage=60
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    log: "logs/reads_profiling/amr/resfinder/{sample}.quantification.log"
    benchmark: "benchmarks/reads_profiling/amr/resfinder/{sample}.quantification.tsv"
    shell:
        """
        bash {input.script} --res {input.res} --mapstat {input.mapstat} --output-prefix {params.prefix} \
            --min-identity {params.min_identity} --min-coverage {params.min_coverage} > {log} 2>&1
        """


rule quantify_virulence_reads:
    input:
        res="data/reads_profiling/virulence/vfdb/{sample}/{sample}.res",
        mapstat="data/reads_profiling/virulence/vfdb/{sample}/{sample}.mapstat",
        script="scripts/quantify_kma_reads.sh"
    output:
        all_refs="data/reads_profiling/virulence/vfdb/{sample}/{sample}.quantification.all.tsv",
        filtered="data/reads_profiling/virulence/vfdb/{sample}/{sample}.quantification.filtered.tsv"
    params:
        prefix="data/reads_profiling/virulence/vfdb/{sample}/{sample}.quantification",
        min_identity=90,
        min_coverage=60
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    log: "logs/reads_profiling/virulence/vfdb/{sample}.quantification.log"
    benchmark: "benchmarks/reads_profiling/virulence/vfdb/{sample}.quantification.tsv"
    shell:
        """
        bash {input.script} --res {input.res} --mapstat {input.mapstat} --output-prefix {params.prefix} \
            --min-identity {params.min_identity} --min-coverage {params.min_coverage} > {log} 2>&1
        """

# Aggregate ResFinder quantifications and annotations.
rule resfinder_summary:
    input:
        quantifications=expand("data/reads_profiling/amr/resfinder/{sample}/{sample}.quantification.filtered.tsv", sample=SAMPLES),
        phenotypes="/archive/extra/daniel.remondini2/bioinf_db/databases/resfinder/2.6.0/resfinder_db/phenotypes.txt"
    output:
        abundance="tables/reads_profiling/resfinder_abundance.tsv",
        metadata="tables/reads_profiling/resfinder_metadata.tsv"
    resources: mem_mb=4000, runtime=30
    container: config["containers"]["r"]
    script: "../scripts/aggregate_resfinder.R"

# Aggregate VFDB quantifications and annotations.
rule vfdb_summary:
    input:
        quantifications=expand("data/reads_profiling/virulence/vfdb/{sample}/{sample}.quantification.filtered.tsv", sample=SAMPLES),
        fasta="/archive/extra/daniel.remondini2/bioinf_db/databases/vfdb/2026-10-26/VFDB_setA_nt.fas",
        annotations="/archive/extra/daniel.remondini2/bioinf_db/databases/vfdb/2026-10-26/VFs.xls"
    output:
        abundance="tables/reads_profiling/vfdb_abundance.tsv",
        metadata="tables/reads_profiling/vfdb_metadata.tsv"
    resources: mem_mb=8000, runtime=30
    container: config["containers"]["r"]
    script: "../scripts/aggregate_vfdb.R"
