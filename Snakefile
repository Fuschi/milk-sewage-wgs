import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_csv(config["samples"], sep="\t", dtype=str)
SAMPLES = samples["sample_id"].tolist()
SAMPLES_BY_BIOME = samples.groupby("biome")["sample_id"].agg(list).to_dict()

# Define workflow targets.
rule all:
    input:
        expand("data/reads_clean/{sample}_{read}_clean.fastq.gz", sample=SAMPLES, read=["R1", "R2", "sing"]),
        expand("tables/read_stats/seqkit_{stage}.tsv", stage=["raw_reads", "trimmed_reads", "cleaned_reads", "cleaned_reads_sing"]),
        expand("data/assembly/single_sample/{sample}/{sample}.contigs.fa", sample=SAMPLES),
        expand("data/assembly/coassembly/{biome}/{biome}.contigs.fa", biome=SAMPLES_BY_BIOME)

include: "rules/cleaning.smk"
include: "rules/stats_reads.smk"
include: "rules/assembly.smk"
include: "rules/coassembly.smk"
