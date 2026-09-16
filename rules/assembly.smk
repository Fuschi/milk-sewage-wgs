# Assemble reads per sample.
rule assembly_single_sample:
    input:
        r1="data/reads_clean/{sample}_R1_clean.fastq.gz",
        r2="data/reads_clean/{sample}_R2_clean.fastq.gz",
        singleton="data/reads_clean/{sample}_sing_clean.fastq.gz",
    output: protected("data/assembly/single_sample/{sample}/{sample}.contigs.fa")
    params: outdir="data/assembly/single_sample/{sample}/megahit"
    log: "logs/assembly/single_sample/{sample}.log"
    benchmark: "benchmarks/assembly/single_sample/{sample}.tsv"
    threads: 8
    resources: mem_mb=32000, runtime=420
    container: config["containers"]["megahit"]
    shell:
        """
        megahit -t {threads} --verbose -1 {input.r1} -2 {input.r2} \
            -r {input.singleton} -o {params.outdir} > {log} 2>&1
        mv {params.outdir}/final.contigs.fa {output}
        """
