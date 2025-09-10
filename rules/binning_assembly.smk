# rules/binnign_assembly.smk
# ─────────────────────────────────────────────
# Binning for assembly conting rules:
# ─────────────────────────────────────────────
#Generic Resources
def default_resources():
    return dict(
        requeue=0,
        trigger=1,
    )
#----------------------------------------------------------------------------------------#
rule map_bbmap_assembly:
    input:
        r1_clean = "snakestream/reads_clean/{sample}_R1_clean.fastq.gz",
        r2_clean = "snakestream/reads_clean/{sample}_R2_clean.fastq.gz",
        sing_clean="snakestream/reads_clean/{sample}_sing_clean.fastq.gz",
        contigs="snakestream/assembly_megahit_genome/{sample_type}/{sample}/{sample}.contigs.fa"
    output:
        outsam="snakestream/map_bbmap/assembly/{sample_type}/{sample}.sam",
    conda: "bbmap"
    log:
        "logs/map_bbmap/assembly/{sample_type}/{sample}.log"
    benchmark:
        "benchmarks/map_bbmap/assembly/{sample_type}/{sample}.txt"
    threads: 20
    resources:
        qos="normal",
        mem_mb=32000,
        time=240,
        **default_resources()
    shell:
            '''
            bbmap.sh \
            in1={input.r1_clean} \
            in2={input.r2_clean} \
            ref={input.contigs} \
            outm={output.outsam} \
            minid=0.90 \
            threads={threads} \
            overwrite=t nodisk=t \
            > {log} 2>&1

 
            '''
#----------------------------------------------------------------------------------------#
rule sort_samtools_assembly:
    input:
        sam="snakestream/map_bbmap/assembly/{sample_type}/{sample}.sam"
    output:
        bam="snakestream/sort_samtools/assembly/{sample_type}/{sample}.bam"
    log:
        "logs/sort_samtools/assembly/{sample_type}/{sample}.log"
    conda:"samtools"
    benchmark:
        "benchmarks/sort_samtools/assembly/{sample_type}/{sample}.txt"
    threads: 4
    resources:
        qos="normal",
        mem_mb=32000,
        time=240,
        **default_resources()
    shell:
        """
        samtools sort -m 20G -@ {threads} -o {output.bam} {input.sam} > {log} 2>&1
        """

#----------------------------------------------------------------------------------------#
rule summarize_depths_assembly:
    input:
        bam="snakestream/sort_samtools/assembly/{sample_type}/{sample}.bam"
    output:
        dep="snakestream/coverage/assembly/{sample_type}/{sample}.txt"
    conda:
        "metabat"
    log:
        "logs/coverage/assembly/{sample_type}/{sample}.log"
    benchmark:
        "benchmarks/coverage/assembly/{sample_type}/{sample}.txt"
    threads: 4
    resources:
        qos="normal",
        mem_mb=16000,
        time=60,
        **default_resources()
    shell:
        """
        jgi_summarize_bam_contig_depths {input.bam} --outputDepth {output.dep} > {log} 2>&1
        """
#----------------------------------------------------------------------------------------#
rule binning_metabat_assembly:
    input:
        contigs="snakestream/assembly_megahit_genome/{sample_type}/{sample}/{sample}.contigs.fa",
        dep="snakestream/coverage/assembly/{sample_type}/{sample}.txt"
    output:
        bin="snakestream/binning_metabat/assembly/{sample_type}/{sample}/{sample}_bin.1.fa",
    params:
        dir_out="snakestream/binning_metabat/assembly/{sample_type}/{sample}",
    conda: "metabat"
    log:
        "logs/binning_metabat/assembly/{sample_type}/{sample}.log"
    benchmark:
        "benchmarks/binning_metabat/assembly/{sample_type}/{sample}.txt"
    threads: 8
    resources:
        qos="normal",
        mem_mb=32000,
        time=240,
        **default_resources()
    shell:
        """
        metabat2 \
        -t {threads} \
        --verbose \
        -i {input.contigs} \
        -a {input.dep} \
        -o {params.dir_out}/{wildcards.sample}_bin \
        -m 1500
        """
 
#----------------------------------------------------------------------------------------#
rule assessment_checkm2_assembly:
    input:
        bin="snakestream/binning_metabat/assembly/{sample_type}/{sample}/{sample}_bin.1.fa"
    output:
        report = "snakestream/checkm2/assembly/{sample_type}/{sample}/{sample}_quality_report.tsv"
    params:
        bins_dir = "snakestream/binning_metabat/assembly/{sample_type}/{sample}",
        dir_out = "snakestream/checkm2/assembly/{sample_type}/{sample}",
        report = "snakestream/checkm2/assembly/{sample_type}/{sample}/quality_report.tsv"
    conda:
        "checkm2"
    log:
       "logs/assessment_checkm2/assembly/{sample_type}/{sample}.out"
    benchmark:
       "benchmarks/assessment_checkm2/asssembly/{sample_type}/{sample}.txt"
    threads: 8
    resources:
       qos="normal",
       mem_mb=32000,
       time=120,
        **default_resources()
    shell:
       """
       checkm2 predict \
               --input {params.bins_dir} \
               --output-directory {params.dir_out} \
               --allmodels \
               --extension fa \
               --force \
               --threads {threads} \
               --database_path databases/CheckM2_database/uniref100.KO.1.dmnd \
               > {log} 2>&1
        mv {params.report} {output.report}
       """
#----------------------------------------------------------------------------------------#

