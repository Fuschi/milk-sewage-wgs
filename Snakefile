# Configuration
import pandas as pd

configfile: "config/config.yaml"

metadata = pd.read_csv(config["samples"], sep="\t", dtype=str)
SAMPLES = metadata["sample_id"].tolist()
SAMPLES_BY_BIOME = metadata.groupby("biome")["sample_id"].apply(list).to_dict()
BIOMES = list(SAMPLES_BY_BIOME)

# -----------------------------------------------------------------------------

# Workflow targets
rule all:
    input:
        # Cleaning
        "tables/kneaddata_read_counts.tsv",

        # Sourmash
        expand("results/sourmash/{sm_group}/{sm_metric}/distance.mat.matrix.png",
               sm_group=["all", "milk", "sewage"], sm_metric=["angular", "jaccard"]),
        expand("results/sourmash/{sm_group}/{sm_metric}/distance.mat.dendro.png",
               sm_group=["all", "milk", "sewage"], sm_metric=["angular", "jaccard"]),

        # Assembly QC
        expand("data/assembly/metaquast/sample/{sample}/report.tsv", sample=SAMPLES),
        expand("data/assembly/metaquast/coassembly/{biome}/report.tsv", biome=BIOMES),

        # Mapping: individual assemblies
        expand("data/contig_mapping/sample/{sample}/{sample}.sorted.bam", sample=SAMPLES),
        expand("data/contig_mapping/sample/{sample}/{sample}.sorted.bam.bai", sample=SAMPLES),
        expand("data/contig_mapping/sample/{sample}/{sample}.coverage.tsv", sample=SAMPLES),

        # Mapping: biome coassemblies
        [f"data/contig_mapping/coassembly/{biome}/{sample}.sorted.bam"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/coassembly/{biome}/{sample}.sorted.bam.bai"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/coassembly/{biome}/{sample}.coverage.tsv"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],

        # Mapping: concatenated assemblies
        [f"data/contig_mapping/multisample/{biome}/{sample}.sorted.bam"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/multisample/{biome}/{sample}.sorted.bam.bai"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/multisample/{biome}/{sample}.coverage.tsv"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],

        # Binning: COMEBin
        expand("data/binning/comebin/sample/{sample}/status.tsv", sample=SAMPLES),
        expand("data/binning/comebin/coassembly/{biome}/comebin_res/comebin_res.tsv", biome=BIOMES),
        expand("data/binning/comebin/coassembly/{biome}/comebin_res/comebin_res_bins", biome=BIOMES),

        # Binning: SemiBin2
        expand("data/binning/semibin2/multisample/{biome}/status.tsv", biome=BIOMES),
        expand("data/binning/semibin2/coassembly/{biome}/output_bins", biome=BIOMES),

        # Binning: VAMB
        expand("data/binning/vamb/multisample/{biome}/status.tsv", biome=BIOMES),
        expand("data/binning/vamb/multisample/{biome}/run", biome=BIOMES)

# -----------------------------------------------------------------------------

# Workflow rules
include: "rules/cleaning.smk"
include: "rules/sourmash.smk"
include: "rules/assembly.smk"
include: "rules/contig_mapping.smk"
include: "rules/binning.smk"
