import pandas as pd

configfile: "config/config.yaml"

metadata = pd.read_csv(config["samples"], sep="\t", dtype=str)
SAMPLES = metadata["sample_id"].tolist()
SAMPLES_BY_BIOME = metadata.groupby("biome")["sample_id"].apply(list).to_dict()
BIOMES = list(SAMPLES_BY_BIOME)

rule all:
    input:
        expand("data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz", sample=SAMPLES),
        expand("data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz", sample=SAMPLES),
        expand("data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz", sample=SAMPLES),
        expand("data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz", sample=SAMPLES),
        "tables/kneaddata_read_counts.tsv",
        expand("data/assembly/sample/{sample}/{sample}.contigs.fa", sample=SAMPLES),
        expand("data/assembly/coassembly/{biome}/{biome}.contigs.fa", biome=BIOMES),
        expand("data/assembly/sample/{sample}/metaquast/report.tsv", sample=SAMPLES),
        expand("data/assembly/coassembly/{biome}/metaquast/report.tsv", biome=BIOMES)

include: "rules/cleaning.smk"
include: "rules/assembly.smk"
