# Assemble each sample independently using paired and unmatched KneadData reads.
# Prefix contig identifiers directly in the MEGAHIT output FASTA.
rule sample_assembly:
    input:
        r1="data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz",
        r2="data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz",
        sing1="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz",
        sing2="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz"
    output:
        protected("data/assembly/sample/{sample}/{sample}.contigs.fa")
    params:
        outdir="data/assembly/sample/{sample}",
        prefix="{sample}",
        memory=30000000000
    threads: 16
    resources:
        mem_mb=32000,
        runtime=480,
        slurm_partition="m1,m2,m3,m4,m5,m6,m7"
    container:
        config["containers"]["megahit"]
    log:
        "logs/assembly/sample/{sample}.log"
    benchmark:
        "benchmarks/assembly/sample/{sample}.tsv"
    shell:
        """
        rm -rf {params.outdir}
        megahit -1 {input.r1} -2 {input.r2} -r {input.sing1},{input.sing2} -t {threads} --memory {params.memory} --verbose \
            -o {params.outdir} --out-prefix {params.prefix} > {log} 2>&1
        sed -i 's/^>/>{wildcards.sample}_/' {output}
        """

# Co-assemble all samples belonging to the same biome using the MEGAHIT meta-large preset.
# Prefix contig identifiers directly in the MEGAHIT output FASTA.
rule coassembly:
    input:
        r1=lambda wc: expand("data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz", sample=SAMPLES_BY_BIOME[wc.biome]),
        r2=lambda wc: expand("data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz", sample=SAMPLES_BY_BIOME[wc.biome]),
        sing1=lambda wc: expand("data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz", sample=SAMPLES_BY_BIOME[wc.biome]),
        sing2=lambda wc: expand("data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz", sample=SAMPLES_BY_BIOME[wc.biome])
    output:
        protected("data/assembly/coassembly/{biome}/{biome}.contigs.fa")
    params:
        r1=lambda wc, input: ",".join(input.r1),
        r2=lambda wc, input: ",".join(input.r2),
        singles=lambda wc, input: ",".join(list(input.sing1) + list(input.sing2)),
        outdir="data/assembly/coassembly/{biome}",
        prefix="{biome}",
        memory=850000000000
    threads: 64
    resources:
        mem_mb=900000,
        runtime=1400,
        slurm_partition="m8",
        slurm_extra="--reservation=prj-trigger --nodelist=mtx30"
    container:
        config["containers"]["megahit"]
    log:
        "logs/assembly/coassembly/{biome}.log"
    benchmark:
        "benchmarks/assembly/coassembly/{biome}.tsv"
    shell:
        """
        if [ -f {params.outdir}/options.json ] && [ -f {params.outdir}/checkpoints.txt ]; then
            echo "===== RESUMING MEGAHIT $(date) =====" >> {log}
            megahit --continue -o {params.outdir} >> {log} 2>&1
        else
            echo "===== STARTING MEGAHIT $(date) =====" >> {log}
            megahit -1 {params.r1} -2 {params.r2} -r {params.singles} -t {threads} --memory {params.memory} --presets meta-large --verbose \
                -o {params.outdir} --out-prefix {params.prefix} >> {log} 2>&1
        fi

        sed -i 's/^>/>{wildcards.biome}_/' {output}
        """

# Assess sample assemblies with MetaQUAST using no-reference settings.
rule sample_contig_qc:
    input:
        "data/assembly/sample/{sample}/{sample}.contigs.fa"
    output:
        "data/assembly/metaquast/sample/{sample}/report.tsv"
    params:
        outdir="data/assembly/metaquast/sample/{sample}"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120
    container:
        config["containers"]["quast"]
    log:
        "logs/assembly/sample/{sample}_metaquast.log"
    benchmark:
        "benchmarks/assembly/sample/{sample}_metaquast.tsv"
    shell:
        """
        metaquast.py {input} -o {params.outdir} -t {threads} \
            --min-contig 200 --max-ref-number 0 --space-efficient > {log} 2>&1
        """


# Assess biome coassemblies with the same no-reference MetaQUAST settings.
rule coassembly_contig_qc:
    input:
        "data/assembly/coassembly/{biome}/{biome}.contigs.fa"
    output:
        "data/assembly/metaquast/coassembly/{biome}/report.tsv"
    params:
        outdir="data/assembly/metaquast/coassembly/{biome}"
    threads: 16
    resources:
        mem_mb=32000,
        runtime=240
    container:
        config["containers"]["quast"]
    log:
        "logs/assembly/coassembly/{biome}_metaquast.log"
    benchmark:
        "benchmarks/assembly/coassembly/{biome}_metaquast.tsv"
    shell:
        """
        metaquast.py {input} -o {params.outdir} -t {threads} \
            --min-contig 200 --max-ref-number 0 --space-efficient > {log} 2>&1
        """
