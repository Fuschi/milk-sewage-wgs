# rules/assembly.smk
# ─────────────────────────────────────────────
# Assembly and Coassembly rules:
# ─────────────────────────────────────────────
#Generic Resources
def default_resources():
    return dict(
        requeue=0,
        trigger=1,
    )
#----------------------------------------------------------------------------------------#
rule assembly_megahit_genome:
    input:
        r1_clean = "snakestream/reads_clean/{sample}_R1_clean.fastq.gz",
        r2_clean = "snakestream/reads_clean/{sample}_R2_clean.fastq.gz",
        sing_clean="snakestream/reads_clean/{sample}_sing_clean.fastq.gz"
    output:
        contigs=protected("snakestream/assembly_megahit_genome/{sample_type}/{sample}/{sample}.contigs.fa"),
    params:
        dir_out="snakestream/assembly_megahit_genome/{sample_type}/{sample}",
        contigs=temp("snakestream/assembly_megahit_genome/{sample_type}/{sample}/final.contigs.fa")
    conda: "megahit"
    log:
        "logs/assembly_megahit_genome/{sample_type}/{sample}.log"
    benchmark:
        "benchmarks/assembly_megahit_genome/{sample_type}/{sample}.txt"
    threads: 8
    resources:
        qos="normal",
        mem_mb=32000,
        time=420,
        **default_resources()
    shell:
        """
        rm -rf {params.dir_out}
        megahit \
        -t {threads} \
        --verbose \
        --min-contig-len 1000 \
        -1 {input.r1_clean} -2 {input.r2_clean} \
        -r {input.sing_clean}\
        -o {params.dir_out} \
           > {log} 2>&1
        mv {params.contigs} {output.contigs}
        """
#--------------------------------------------------------------------------------------#
rule reads_concatenation:
    input:
        r1_clean=lambda wildcards: expand("snakestream/reads_clean/{sample}_R1_clean.fastq.gz",sample=BIOME_TO_SAMPLE[wildcards.sample_type]),
        r2_clean=lambda wildcards: expand("snakestream/reads_clean/{sample}_R2_clean.fastq.gz",sample=BIOME_TO_SAMPLE[wildcards.sample_type]),
        sing_clean=lambda wildcards: expand("snakestream/reads_clean/{sample}_sing_clean.fastq.gz",sample=BIOME_TO_SAMPLE[wildcards.sample_type]),
    output:
        conc_R1="snakestream/coassembly_megahit_genome/conc_reads/{sample_type}/all_R1.fastq.gz",
        conc_R2="snakestream/coassembly_megahit_genome/conc_reads/{sample_type}/all_R2.fastq.gz",
        conc_sing="snakestream/coassembly_megahit_genome/conc_reads/{sample_type}/all_sing.fastq.gz",
    conda: "megahit"
    log:
        out="logs/reads_concatenation/{sample_type}.out",
        err="logs/reads_concatenation/{sample_type}.err"
    benchmark:
        "benchmarks/reads_concatenation/{sample_type}.txt"
    threads: 32
    resources:
        qos="normal",
        mem_mb=500000,
        time=1430,
        **default_resources()
    shell:
        """
        cat {input.r1_clean} | gzip -c > {output.conc_R1}

        cat {input.r2_clean} | gzip -c > {output.conc_R2}

        cat {input.sing_clean} | gzip -c > {output.conc_sing}
        """
##NOT WORIKING
#--------------------------------------------------------------------------------------#
rule coassembly_megahit_genome:
    input:
        conc_R1="snakestream/coassembly_megahit_genome/conc_reads/{sample_type}/all_R1.fastq.gz",
        conc_R2="snakestream/coassembly_megahit_genome/conc_reads/{sample_type}/all_R2.fastq.gz",
        conc_sing="snakestream/coassembly_megahit_genome/conc_reads/{sample_type}/all_sing.fastq.gz",
    output:
        contigs=protected("snakestream/coassembly_megahit_genome/{sample_type}/coassembly.final.contigs.fa")
    params:
        dir_out="snakestream/coassembly_megahit_genome/{sample_type}",
        contigs=temp("snakestream/coassembly_megahit_genome/{sample_type}/final.contigs.fa")
    conda: "megahit"
    log:
        out="logs/coassembly_megahit_genome/{sample_type}.out",
        err="logs/coassembly_megahit_genome/{sample_type}.err"
    benchmark:
        "benchmarks/coassembly_megahit_genome/{sample_type}.txt"
    threads: 64
    resources:
        qos="long",
        mem_mb=1000000,
        time=1430,
        requeue=1,
        trigger=1
    shell:
        """
        megahit \
            -t {threads} \
            --verbose \
            -1 {input.conc_R1} -2 {input.conc_R2} \
            -r {input.conc_sing} \
            -o {params.dir_out}
        mv {params.contigs} {output.contigs}
        """
#--------------------------------------------------------------------------------------#
