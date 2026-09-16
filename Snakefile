import pandas as pd

configfile: "config/config.yaml"

SAMPLES = pd.read_csv(config["samples"], sep="\t", dtype=str)["sample_id"].tolist()

# Define workflow targets.
rule all:
    input:
        expand("data/reads_clean/{sample}_{read}_clean.fastq.gz", sample=SAMPLES, read=["R1", "R2", "sing"]),
        expand("tables/read_stats/seqkit_{stage}.tsv", stage=["raw_reads", "trimmed_reads", "cleaned_reads", "cleaned_reads_sing"])

include: "rules/cleaning.smk"
include: "rules/stats_reads.smk"
