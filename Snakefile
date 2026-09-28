import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_csv(config["samples"], sep="\t", dtype=str, keep_default_na=False)
if "sample_id" not in samples.columns or samples.empty:
    raise ValueError("The sample sheet must have a sample_id column and at least one sample.")
if samples["sample_id"].duplicated().any():
    raise ValueError("Duplicate sample_id values in the sample sheet.")
if not samples["sample_id"].str.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*").all():
    raise ValueError("Sample IDs must contain only letters, digits, underscores or hyphens.")
SAMPLES = samples["sample_id"].tolist()

rule all:
    input:
        expand("data/kneaddata/{sample}", sample=SAMPLES),
        expand("tables/read_stats/kneaddata/{sample}.tsv", sample=SAMPLES),

include: "rules/cleaning.smk"
