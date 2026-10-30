SOURMASH_GROUPS = {
    "all": SAMPLES,
    "milk": SAMPLES_BY_BIOME["milk"],
    "sewage": SAMPLES_BY_BIOME["sewage"]
}

SOURMASH_METRICS = {
    "angular": "",
    "jaccard": "--ignore-abundance"
}


# Build one abundance-aware Sourmash sketch per sample from all cleaned reads.
rule sourmash_sketch:
    input:
        r1="data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz",
        r2="data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz",
        sing1="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz",
        sing2="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz"
    output: "data/sourmash_sketches/{sample}.zip"
    threads: 1
    resources: mem_mb=4000, runtime=60
    container: config["containers"]["sourmash"]
    log: "logs/sourmash/sketch/{sample}.log"
    benchmark: "benchmarks/sourmash/sketch/{sample}.tsv"
    shell:
        """
        sourmash sketch dna -p k=31,scaled=1000,abund --merge {wildcards.sample} \
            -o {output} {input.r1} {input.r2} {input.sing1} {input.sing2} > {log} 2>&1
        """


# Compare all, milk, or sewage samples using abundance-aware angular or Jaccard distance.
rule sourmash_distance:
    input: lambda wc: expand("data/sourmash_sketches/{sample}.zip", sample=SOURMASH_GROUPS[wc.sm_group])
    output:
        csv="results/sourmash/{sm_group}_{sm_metric}_distance.csv",
        matrix="results/sourmash/{sm_group}/{sm_metric}/distance.mat",
        labels="results/sourmash/{sm_group}/{sm_metric}/labels.csv"
    params: metric=lambda wc: SOURMASH_METRICS[wc.sm_metric]
    threads: 1
    resources: mem_mb=8000, runtime=60
    container: config["containers"]["sourmash"]
    log: "logs/sourmash/{sm_group}_{sm_metric}.log"
    benchmark: "benchmarks/sourmash/{sm_group}_{sm_metric}.tsv"
    shell:
        """
        sourmash compare {input} {params.metric} --distance-matrix \
            --csv {output.csv} --output {output.matrix} --labels-to {output.labels} > {log} 2>&1
        """

# Cluster each Sourmash distance matrix and generate heatmap, dendrogram, and histogram.
rule sourmash_plot:
    input:
        matrix="results/sourmash/{sm_group}/{sm_metric}/distance.mat",
        labels="results/sourmash/{sm_group}/{sm_metric}/labels.csv"
    output:
        matrix_png="results/sourmash/{sm_group}/{sm_metric}/distance.mat.matrix.png",
        dendrogram_png="results/sourmash/{sm_group}/{sm_metric}/distance.mat.dendro.png",
        hist_png="results/sourmash/{sm_group}/{sm_metric}/distance.mat.hist.png"
    threads: 1
    resources: mem_mb=4000, runtime=30
    container: config["containers"]["sourmash"]
    log: "logs/sourmash/{sm_group}_{sm_metric}_plot.log"
    benchmark: "benchmarks/sourmash/{sm_group}_{sm_metric}_plot.tsv"
    shell:
        """
        ROOT=$PWD
        cd results/sourmash/{wildcards.sm_group}/{wildcards.sm_metric}
        sourmash plot distance.mat --labels-from labels.csv > "$ROOT/{log}" 2>&1
        """

