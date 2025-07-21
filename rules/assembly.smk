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
        time=240,
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
##NOT WORIKING
rule coassembly_megahit_genome:
    input:
        R1="snakestream/reads_clean/{sample_type}/{sample}_R1.clean.fastq.gz",
        R2="snakestream/reads_clean/{sample_type}/{sample}_R2.clean.fastq.gz",
    output:
        contigs=protected("snakestream/coassembly_megahit_genome/{sample_type}/coassembly.final.contigs.fa")
    params:
        dir_out="coassembly_megahit_genome/{sample_type}",
        tmp_dir=temp("coassembly_megahit_genome/tmp/{sample_type}"),
        temp_R1=temp("coassembly_megahit_genome/tmp/{sample_type}/all_R1.fastq.gz"),
        temp_R2=temp("coassembly_megahit_genome/tmp/{sample_type}/all_R2.fastq.gz"),
        contigs=temp("coassembly_megahit_genome/{sample_type}/final.contigs.fa")
    conda: "megahit"
    log:
        out="logs/coassembly_megahit_genome/{sample_type}.out",
        err="logs/coassembly_megahit_genome/{sample_type}.err"
    benchmark:
        "benchmarks/coassembly_megahit_genome/{sample_type}.txt"
    threads: 64
    resources:
        qos="normal",
        mem_mb=192000,
        time=1000,
        **default_resources()
    shell:
        """
        mkdir -p {params.tmp_dir}
        rm -rf {params.dir_out}
        cat {input.R1} > {params.temp_R1}

        cat {input.R2} > {params.temp_R2}

        megahit \
            -t {threads} \
            --verbose \
            --min-contig-len 1000 \
            -1 {params.temp_R1} -2 {params.temp_R2} \
            -o {params.dir_out}
        rm -rf {params.tmp_dir}
        mv {params.contigs} {output.contigs}
        """
#--------------------------------------------------------------------------------------#
