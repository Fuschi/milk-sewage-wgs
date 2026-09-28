# Bos taurus Bowtie2 index used by KneadData.
HOST_INDEX = config["references"]["bostaurus_index"]

HOST_INDEX_FILES = [
    f"{HOST_INDEX}.{suffix}.bt2"
    for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")
]


# Clean paired-end reads with KneadData:
# - FastQC before cleaning
# - adapter/quality trimming with Trimmomatic
# - Bos taurus host removal with Bowtie2
# - FastQC after cleaning
#
# KneadData also manages reads that become unmatched during trimming/filtering.
rule kneaddata_clean:
    input:
        r1="data/reads_raw/{sample}_R1.fastq.gz",
        r2="data/reads_raw/{sample}_R2.fastq.gz",
        index=HOST_INDEX_FILES
    output:
        reads=directory("data/kneaddata/{sample}")
    params:
        index=HOST_INDEX,
        prefix="{sample}_kneaddata"
    threads: 8
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
            --output {output.reads} \
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
            > {log} 2>&1
        """


# Generate a single read-count summary for all samples.
rule kneaddata_read_counts:
    input:
        logs=expand("logs/kneaddata/{sample}.log", sample=SAMPLES)
    output:
        "tables/read_stats/kneaddata_read_counts.tsv"
    resources:
        mem_mb=2000,
        runtime=10
    container:
        config["containers"]["kneaddata"]
    log:
        "logs/kneaddata/read_counts.log"
    benchmark:
        "benchmarks/kneaddata/read_counts.tsv"
    shell:
        """
        kneaddata_read_count_table \
            --input {input.logs} \
            --output {output} \
            > {log} 2>&1
        """