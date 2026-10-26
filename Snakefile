# Configuration
import pandas as pd

configfile: "config/config.yaml"

metadata = pd.read_csv(config["samples"], sep="\t", dtype=str)
SAMPLES = metadata["sample_id"].tolist()
SAMPLES_BY_BIOME = metadata.groupby("biome")["sample_id"].apply(list).to_dict()
BIOMES = list(SAMPLES_BY_BIOME)
BIOME_BY_SAMPLE = metadata.set_index("sample_id")["biome"].to_dict()

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

        # Reads profiling: AMR
        expand("data/reads_profiling/amr/resfinder/{sample}/{sample}.res", sample=SAMPLES),
        expand("data/reads_profiling/amr/resfinder/{sample}/{sample}.mapstat", sample=SAMPLES),
        # Reads profiling: Virulence
        expand("data/reads_profiling/amr/resfinder/{sample}/{sample}.quantification.{table}.tsv", sample=SAMPLES, table=["all", "filtered"]),
        expand("data/reads_profiling/virulence/vfdb/{sample}/{sample}.quantification.{table}.tsv", sample=SAMPLES, table=["all", "filtered"]),
        # Read profiling: summaries
        expand("tables/reads_profiling/{database}_{table}.tsv", database=["resfinder", "vfdb"], table=["abundance", "metadata"]),


        # Assembly QC
        expand("data/assembly/metaquast/sample/{sample}/report.tsv", sample=SAMPLES),
        expand("data/assembly/metaquast/coassembly/{biome}/report.tsv", biome=BIOMES),

        # Assembly QC summary
        "tables/metaquast_summary.tsv",

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
        # COMEBin filenames
        expand("data/binning/comebin/sample/{sample}/bins_renamed.done", sample=SAMPLES),
        expand("data/binning/comebin/coassembly/{biome}/bins_renamed.done", biome=BIOMES),

        # Binning: SemiBin2
        expand("data/binning/semibin2/multisample/{biome}/status.tsv", biome=BIOMES),
        expand("data/binning/semibin2/coassembly/{biome}/output_bins", biome=BIOMES),
        # SemiBin2 headers
        expand("data/binning/semibin2/{mode}/{biome}/headers_restored.done", mode=["multisample", "coassembly"], biome=BIOMES),

        # Binning: VAMB
        expand("data/binning/vamb/multisample/{biome}/status.tsv", biome=BIOMES),
        expand("data/binning/vamb/multisample/{biome}/run", biome=BIOMES),
        # VAMB filenames
        expand("data/binning/vamb/multisample/{biome}/bins_renamed.done", biome=BIOMES),

        # Refinement bin inputs
        expand("data/refinement/dastool/sample/{sample}/inputs/{tool}.tsv", sample=SAMPLES, tool=["comebin", "semibin2", "vamb"]),
        expand("data/refinement/dastool/coassembly/{biome}/inputs/{tool}.tsv", biome=BIOMES, tool=["comebin", "semibin2"]),

        # Refinement: sample bins
        expand("data/refinement/dastool/sample/{sample}/{sample}_DASTool_summary.tsv", sample=SAMPLES),
        expand("data/refinement/dastool/sample/{sample}/{sample}_DASTool_contig2bin.tsv", sample=SAMPLES),
        expand("data/refinement/dastool/sample/{sample}/{sample}_DASTool_bins", sample=SAMPLES),
        expand("data/refinement/dastool/sample/{sample}/status.tsv", sample=SAMPLES),

        # Refinement: coassembly bins
        expand("data/refinement/dastool/coassembly/{biome}/{biome}_DASTool_summary.tsv", biome=BIOMES),
        expand("data/refinement/dastool/coassembly/{biome}/{biome}_DASTool_contig2bin.tsv", biome=BIOMES),
        expand("data/refinement/dastool/coassembly/{biome}/{biome}_DASTool_bins", biome=BIOMES),

        # Bin quality: samples
        expand("data/bin_quality/checkm2/sample/{sample}/quality_report.tsv", sample=SAMPLES),
        expand("data/bin_quality/gunc/sample/{sample}/GUNC.progenomes_3.maxCSS_level.tsv", sample=SAMPLES),

        # Bin quality: coassemblies
        expand("data/bin_quality/checkm2/coassembly/{biome}/quality_report.tsv", biome=BIOMES),
        expand("data/bin_quality/gunc/coassembly/{biome}/GUNC.progenomes_3.maxCSS_level.tsv", biome=BIOMES),

# -----------------------------------------------------------------------------

# Workflow rules
include: "rules/cleaning.smk"
include: "rules/sourmash.smk"
include: "rules/reads_profiling.smk"
include: "rules/assembly.smk"
include: "rules/contig_mapping.smk"
include: "rules/binning.smk"
include: "rules/refinement.smk"
include: "rules/bin_quality.smk"
