import pandas as pd

configfile: "config/config.yaml"

metadata = pd.read_csv(config["samples"], sep="\t", dtype=str)
SAMPLES = metadata["sample_id"].tolist()
SAMPLES_BY_BIOME = metadata.groupby("biome")["sample_id"].apply(list).to_dict()
BIOMES = list(SAMPLES_BY_BIOME)

rule all:
    input:
        # cleaning
        "tables/kneaddata_read_counts.tsv",

        # sourmash
        expand("results/sourmash/{sm_group}/{sm_metric}/distance.mat.matrix.png",
               sm_group=["all", "milk", "sewage"], sm_metric=["angular", "jaccard"]),
        expand("results/sourmash/{sm_group}/{sm_metric}/distance.mat.dendro.png",
               sm_group=["all", "milk", "sewage"], sm_metric=["angular", "jaccard"]),

        # assembly
        expand("data/assembly/metaquast/sample/{sample}/report.tsv", sample=SAMPLES),
        expand("data/assembly/metaquast/coassembly/{biome}/report.tsv", biome=BIOMES),

        # mapping - individual assemblies
        expand("data/contig_mapping/sample/{sample}/{sample}.sorted.bam", sample=SAMPLES),
        expand("data/contig_mapping/sample/{sample}/{sample}.sorted.bam.bai", sample=SAMPLES),
        expand("data/contig_mapping/sample/{sample}/{sample}.coverage.tsv", sample=SAMPLES),

        # mapping - biome coassemblies
        [f"data/contig_mapping/coassembly/{biome}/{sample}.sorted.bam"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/coassembly/{biome}/{sample}.sorted.bam.bai"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/coassembly/{biome}/{sample}.coverage.tsv"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],

        # mapping - concatenated individual assemblies
        [f"data/contig_mapping/multisample/{biome}/{sample}.sorted.bam"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/multisample/{biome}/{sample}.sorted.bam.bai"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/multisample/{biome}/{sample}.coverage.tsv"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],

        # binning - COMEBin
        expand("data/binning/comebin/sample/{sample}/status.tsv", sample=SAMPLES),
        # binning - SemiBin2
        expand("data/binning/semibin2/multisample/{biome}/status.tsv", biome=BIOMES),
        # binning - VAMB
        expand("data/binning/vamb/multisample/{biome}/status.tsv", biome=BIOMES),
        expand("data/binning/vamb/multisample/{biome}/run", biome=BIOMES)

include: "rules/cleaning.smk"
include: "rules/sourmash.smk"
include: "rules/assembly.smk"
include: "rules/contig_mapping.smk"
include: "rules/binning.smk"
