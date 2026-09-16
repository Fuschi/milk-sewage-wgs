import pandas as pd

configfile: "config/config.yaml"

SAMPLES = pd.read_csv(config["samples"], sep="\t", dtype=str)["sample_id"].tolist()

# Define workflow targets.
rule all:
    input: expand("data/reads_clean/{sample}_{read}_clean.fastq.gz", sample=SAMPLES, read=["R1", "R2", "sing"])

include: "rules/cleaning.smk"
