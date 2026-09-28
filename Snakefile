import pandas as pd

configfile: "config/config.yaml"

metadata = pd.read_csv(config["samples"], sep="\t", dtype=str)
SAMPLES = metadata["sample_id"].tolist()

rule all:
    input:
        expand("data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz", sample=SAMPLES),
        expand("data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz", sample=SAMPLES),
        expand("data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz", sample=SAMPLES),
        expand("data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz", sample=SAMPLES),
        "tables/kneaddata_read_counts.tsv"

include: "rules/cleaning.smk"
