# rules/binning_coassembly.smk
# ─────────────────────────────────────────────
# Binning for Coassembly contigs rules:
# ─────────────────────────────────────────────
#Generic Resources
def default_resources():
    return dict(
        requeue=0,
        trigger=1,
    )
#----------------------------------------------------------------------------------------#
rule map_bbmap_coassembly:
    input:
        conc_R1="snakestream/coassembly_megahit_genome/conc_reads/{sample_type}/all_R1.fastq.gz",
        conc_R2="snakestream/coassembly_megahit_genome/conc_reads/{sample_type}/all_R2.fastq.gz",
        conc_sing="snakestream/coassembly_megahit_genome/conc_reads/{sample_type}/all_sing.fastq.gz",
        contigs="snakestream/coassembly_megahit_genome/{sample_type}/coassembly.final.contigs.fa"
    output:
        outsam="snakestream/map_bbmap/coassembly/{sample_type}.sam",
    conda: "bbmap"
    log:
        "logs/map_bbmap/coassembly/{sample_type}.log"
    benchmark:
        "benchmarks/map_bbmap/coassembly/{sample_type}.txt"
    threads: 20
    resources:
        qos="normal",
        mem_mb=32000,
        time=240,
        **default_resources()
    shell:
            '''
            bbmap.sh \
            in1={input.conc_R1} \
            in2={input.conc_R2} \
            ref={input.contigs} \
            outm={output.outsam} \
            minid=0.90 \
            threads={threads} \
            overwrite=t nodisk=t \
            > {log} 2>&1
            '''
#----------------------------------------------------------------------------------------#
rule sort_samtools_coassembly:
    input:
        sam="snakestream/map_bbmap/coassembly/{sample_type}.sam"
    output:
        bam="snakestream/sort_samtools/coassembly/{sample_type}.bam"
    log:
        "logs/sort_samtools/coassembly/{sample_type}.log"
    conda:"samtools"
    benchmark:
        "benchmarks/sort_samtools/coassembly/{sample_type}.txt"
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
rule summarize_depths_coassembly:
    input:
        bam="snakestream/sort_samtools/assembly/{sample_type}.bam"
    output:
        dep="snakestream/coverage/coassembly/{sample_type}.txt"
    conda:
        "metabat"
    log:
        "logs/coverage/coassembly/{sample_type}.log"
    benchmark:
        "benchmarks/coverage/coassembly/{sample_type}.txt"
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
rule binning_metabat_coassembly:
    input:
        contigs="snakestream/coassembly_megahit_genome/{sample_type}/coassembly.final.contigs.fa",
        dep="snakestream/coverage/coassembly/{sample_type}.txt"
    output:
        bin="snakestream/binning_metabat/coassembly/{sample_type}/{sample_type}_bin.1.fa",
    params:
        dir_out="snakestream/binning_metabat/coassembly/{sample_type}",
    conda: "metabat"
    log:
        "logs/binning_metabat/coassembly/{sample_type}.log"
    benchmark:
        "benchmarks/binning_metabat/coassembly/{sample_type}.txt"
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
        -o {params.dir_out}/{wildcards.sample_type}_bin \
        -m 1500
        > {log} 2>&1
        """

#----------------------------------------------------------------------------------------#
rule assessment_checkm2_coassembly:
    input:
        bin="snakestream/binning_metabat/coassembly/{sample_type}/{sample_type}_bin.1.fa"
    output:
        report = "snakestream/checkm2/coassembly/{sample_type}/{sample_type}_quality_report.tsv"
    params:
        bins_dir = "snakestream/binning_metabat/coassembly/{sample_type}",
        dir_out = "snakestream/checkm2/coassembly/{sample_type}/",
        report = "snakestream/checkm2/coassembly/{sample_type}/quality_report.tsv"
    conda:
        "checkm2"
    log:
       "logs/assessment_checkm2/coassembly/{sample_type}.log"
    benchmark:
       "benchmarks/assessment_checkm2/coasssembly/{sample_type}.txt"
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
