import pandas as pd

configfile: "config/config.yaml"

metadata = pd.read_csv(config["samples"], sep="\t", dtype=str)
SAMPLES = metadata["sample_id"].tolist()
SAMPLES_BY_BIOME = metadata.groupby("biome")["sample_id"].apply(list).to_dict()
BIOMES = list(SAMPLES_BY_BIOME)
BIOME_BY_SAMPLE = metadata.set_index("sample_id")["biome"].to_dict()

rule all:
    input:
        # Cleaning
        "results/kneaddata_read_counts.tsv",

        # Sourmash
        expand("results/sourmash/{sm_group}/{sm_metric}/distance.mat.matrix.png",
               sm_group=["all", "milk", "sewage"], sm_metric=["angular", "jaccard"]),
        expand("results/sourmash/{sm_group}/{sm_metric}/distance.mat.dendro.png",
               sm_group=["all", "milk", "sewage"], sm_metric=["angular", "jaccard"]),

        # Read profiling
        expand("data/resfinder_reads/{sample}/{sample}.res", sample=SAMPLES),
        expand("data/resfinder_reads/{sample}/{sample}.mapstat", sample=SAMPLES),
        expand("data/resfinder_reads/{sample}/{sample}.quantification.{table}.tsv", sample=SAMPLES, table=["all", "filtered"]),
        expand("data/vfdb_reads/{sample}/{sample}.quantification.{table}.tsv", sample=SAMPLES, table=["all", "filtered"]),
        expand("results/{database}_reads_abundance.tsv", database=["resfinder", "vfdb"]),
        expand("results/{database}_metadata.tsv", database=["resfinder", "vfdb"]),

        # Assembly quality
        expand("data/metaquast/sample/{sample}/report.tsv", sample=SAMPLES),
        expand("data/metaquast/coassembly/{biome}/report.tsv", biome=BIOMES),

        "results/metaquast_summary.tsv",

        # Mapping and coverage
        expand("data/contig_mapping/sample/{sample}/{sample}.sorted.bam", sample=SAMPLES),
        expand("data/contig_mapping/sample/{sample}/{sample}.sorted.bam.bai", sample=SAMPLES),
        expand("data/contig_mapping/sample/{sample}/{sample}.coverage.tsv", sample=SAMPLES),

        [f"data/contig_mapping/coassembly/{biome}/{sample}.sorted.bam"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/coassembly/{biome}/{sample}.sorted.bam.bai"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/coassembly/{biome}/{sample}.coverage.tsv"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],

        [f"data/contig_mapping/multisample/{biome}/{sample}.sorted.bam"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/multisample/{biome}/{sample}.sorted.bam.bai"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],
        [f"data/contig_mapping/multisample/{biome}/{sample}.coverage.tsv"
         for biome, samples in SAMPLES_BY_BIOME.items() for sample in samples],

        # Binning
        expand("data/comebin/sample/{sample}/status.tsv", sample=SAMPLES),
        expand("data/comebin/coassembly/{biome}/comebin_res/comebin_res.tsv", biome=BIOMES),
        expand("data/comebin/coassembly/{biome}/comebin_res/comebin_res_bins", biome=BIOMES),
        expand("data/comebin/sample/{sample}/bins_renamed.done", sample=SAMPLES),
        expand("data/comebin/coassembly/{biome}/bins_renamed.done", biome=BIOMES),

        expand("data/semibin2/multisample/{biome}/status.tsv", biome=BIOMES),
        expand("data/semibin2/coassembly/{biome}/output_bins", biome=BIOMES),
        expand("data/semibin2/{mode}/{biome}/headers_restored.done", mode=["multisample", "coassembly"], biome=BIOMES),

        expand("data/vamb/multisample/{biome}/status.tsv", biome=BIOMES),
        expand("data/vamb/multisample/{biome}/run", biome=BIOMES),
        expand("data/vamb/multisample/{biome}/bins_renamed.done", biome=BIOMES),

        # Bin refinement and quality
        expand("data/dastool/sample/{sample}/inputs/{tool}.tsv", sample=SAMPLES, tool=["comebin", "semibin2", "vamb"]),
        expand("data/dastool/coassembly/{biome}/inputs/{tool}.tsv", biome=BIOMES, tool=["comebin", "semibin2"]),

        expand("data/dastool/sample/{sample}/{sample}_DASTool_summary.tsv", sample=SAMPLES),
        expand("data/dastool/sample/{sample}/{sample}_DASTool_contig2bin.tsv", sample=SAMPLES),
        expand("data/dastool/sample/{sample}/{sample}_DASTool_bins", sample=SAMPLES),
        expand("data/dastool/sample/{sample}/status.tsv", sample=SAMPLES),

        expand("data/dastool/coassembly/{biome}/{biome}_DASTool_summary.tsv", biome=BIOMES),
        expand("data/dastool/coassembly/{biome}/{biome}_DASTool_contig2bin.tsv", biome=BIOMES),
        expand("data/dastool/coassembly/{biome}/{biome}_DASTool_bins", biome=BIOMES),

        expand("data/checkm2/sample/{sample}/quality_report.tsv", sample=SAMPLES),
        expand("data/gunc/sample/{sample}/GUNC.progenomes_3.maxCSS_level.tsv", sample=SAMPLES),

        expand("data/checkm2/coassembly/{biome}/quality_report.tsv", biome=BIOMES),
        expand("data/gunc/coassembly/{biome}/GUNC.progenomes_3.maxCSS_level.tsv", biome=BIOMES),

        "results/bin_quality_summary.tsv",
        "data/drep_inputs/genomes.txt",
        "data/drep_inputs/genome_info.csv",

        # Dereplication, taxonomy, and MAG profiling
        "data/drep",

        "data/gtdbtk",
        "results/gtdbtk_summary.tsv",
        "results/ncbi_taxonomy.tsv",

        expand("data/coverm_mags/{sample}.tsv", sample=SAMPLES),
        expand("results/abricate_{database}.tsv", database=["resfinder_2.6.0", "vfdb_2026-10-26"]),
        "results/mags_abundance.tsv"

include: "rules/cleaning.smk"
include: "rules/sourmash.smk"
include: "rules/reads_profiling.smk"
include: "rules/assembly.smk"
include: "rules/contig_mapping.smk"
include: "rules/binning.smk"
include: "rules/refinement.smk"
include: "rules/bin_quality.smk"
include: "rules/dereplication.smk"
include: "rules/taxonomy_mag.smk"
include: "rules/mags_profiling.smk"
