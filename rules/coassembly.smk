# Assemble reads by biome.
rule coassembly:
    input:
        r1=lambda wc: expand("data/reads_clean/{sample}_R1_clean.fastq.gz", sample=SAMPLES_BY_BIOME[wc.biome]),
        r2=lambda wc: expand("data/reads_clean/{sample}_R2_clean.fastq.gz", sample=SAMPLES_BY_BIOME[wc.biome]),
        singleton=lambda wc: expand("data/reads_clean/{sample}_sing_clean.fastq.gz", sample=SAMPLES_BY_BIOME[wc.biome]),
    output: protected("data/assembly/coassembly/{biome}/{biome}.contigs.fa")
    params:
        r1=lambda wc, input: ",".join(input.r1),
        r2=lambda wc, input: ",".join(input.r2),
        singleton=lambda wc, input: ",".join(input.singleton),
        outdir="data/assembly/coassembly/{biome}/megahit",
        memory=850000000000,
    log: "logs/assembly/coassembly/{biome}.log"
    benchmark: "benchmarks/assembly/coassembly/{biome}.tsv"
    threads: 64
    resources: cpus_per_task=64, mem_mb=900000, runtime=4320, slurm_qos="long", slurm_partition="m8"
    container: config["containers"]["megahit"]
    shell:
        """
        megahit -t {threads} --memory {params.memory} --verbose \
            -1 {params.r1} -2 {params.r2} -r {params.singleton} -o {params.outdir} > {log} 2>&1
        mv {params.outdir}/final.contigs.fa {output}
        """
