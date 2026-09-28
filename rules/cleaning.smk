HOST_INDEX = config["references"]["bostaurus_index"]
HOST_INDEX_FILES = [f"{HOST_INDEX}.{suffix}.bt2" for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")]

rule kneaddata_clean:
    input:
        r1="data/reads_raw/{sample}_R1.fastq.gz",
        r2="data/reads_raw/{sample}_R2.fastq.gz",
        index=HOST_INDEX_FILES
    output:
        r1=temp("data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq"),
        r2=temp("data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq"),
        sing1=temp("data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq"),
        sing2=temp("data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq")
    params:
        outdir="data/kneaddata/{sample}",
        index=HOST_INDEX,
        prefix="{sample}_kneaddata"
    threads: 16
    resources:
        mem_mb=32000,
        runtime=240
    container:
        config["containers"]["kneaddata"]
    log:
        "logs/kneaddata/{sample}.log"
    benchmark:
        "benchmarks/kneaddata/{sample}.tsv"
    shell:
        """
        kneaddata \
            --input1 {input.r1} \
            --input2 {input.r2} \
            --output {params.outdir} \
            --output-prefix {params.prefix} \
            --reference-db {params.index} \
            --threads {threads} \
            --processes 1 \
            --max-memory 4G \
            --sequencer-source NexteraPE \
            --quality-scores phred33 \
            --bowtie2-options="--very-sensitive-local" \
            --decontaminate-pairs strict \
            --reorder \
            --bypass-trf \
            --remove-intermediate-output \
            --run-fastqc-start \
            --run-fastqc-end \
            --log {log}
        """

rule compress_kneaddata_reads:
    input:
        r1="data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq",
        r2="data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq",
        sing1="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq",
        sing2="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq"
    output:
        r1="data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz",
        r2="data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz",
        sing1="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz",
        sing2="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz"
    threads: 8
    resources:
        mem_mb=4000,
        runtime=60
    benchmark:
        "benchmarks/kneaddata/compress_{sample}.tsv"
    shell:
        """
        pigz -p {threads} -c {input.r1} > {output.r1}
        pigz -p {threads} -c {input.r2} > {output.r2}
        pigz -p {threads} -c {input.sing1} > {output.sing1}
        pigz -p {threads} -c {input.sing2} > {output.sing2}
        """

rule kneaddata_read_counts:
    input:
        r1=expand("data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz", sample=SAMPLES),
        r2=expand("data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz", sample=SAMPLES),
        sing1=expand("data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz", sample=SAMPLES),
        sing2=expand("data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz", sample=SAMPLES),
        logs=expand("logs/kneaddata/{sample}.log", sample=SAMPLES)
    output:
        "tables/kneaddata_read_counts.tsv"
    resources:
        mem_mb=2000,
        runtime=10
    container:
        config["containers"]["kneaddata"]
    benchmark:
        "benchmarks/kneaddata/read_counts.tsv"
    shell:
        """
        kneaddata_read_count_table --input logs/kneaddata --output {output}
        """

