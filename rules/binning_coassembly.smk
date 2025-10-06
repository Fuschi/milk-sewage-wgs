# rules/binning_coassembly.smk
# ─────────────────────────────────────────────
import glob
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
    threads: 64
    resources:
        qos="normal",
        mem_mb=500000,
        time=1430,
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
    threads: 32
    resources:
        qos="normal",
        mem_mb=700000,
        time=1430,
        **default_resources()
    shell:
        """
        samtools view -@ {threads} -bS {input.sam} | samtools sort -m 20G -@ {threads} -o {output.bam} > {log} 2>&1
        """

#----------------------------------------------------------------------------------------#
rule summarize_depths_coassembly:
    input:
        bam="snakestream/sort_samtools/coassembly/{sample_type}.bam"
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
        bin="snakestream/binning_metabat/coassembly/{sample_type}/{sample_type}_bin.BinInfo.txt",
    params:
        dir_out="snakestream/binning_metabat/coassembly/{sample_type}",
    conda: "metabat"
    log:
        "logs/binning_metabat/coassembly/{sample_type}.log"
    benchmark:
        "benchmarks/binning_metabat/coassembly/{sample_type}.txt"
    threads: 32
    resources:
        qos="normal",
        mem_mb=32000,
        time=1430,
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
        bin="snakestream/binning_metabat/coassembly/{sample_type}/{sample_type}_bin.BinInfo.txt"
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
rule merge_checkm2_reports_coassembly:
    input:
        report = "snakestream/checkm2/coassembly/{sample_type}/{sample_type}_quality_report.tsv"
    output:
        tables="snakestream/tables/coassembly/{sample_type}/checkm2_reports_for_dRep.csv"
    resources:
        qos="normal",
        mem_mb=32000,
        time=120
    conda:
        "r-tidyverse"
    script:
        "scripts/merge_checkm2_reports_coassembly.R"


#----------------------------------------------------------------------------------------#
rule dereplicate_genomes_assembly:
    input:
          bins="snakestream/binning_metabat/coassembly/{sample_type}/{sample_type}_bin.BinInfo.txt",
          checkm2_table="snakestream/tables/coassembly/{sample_type}/checkm2_reports_for_dRep.csv"
    output:
        "snakestream/dereplicated_genome/coassembly/{sample_type}/figures/Winning_genomes.pdf"
    params:
        dir_out=lambda wildcards: f"snakestream/dereplicated_genome/coassembly/{wildcards.sample_type}/",
        bins=lambda wildcards: sorted(glob(f"snakestream/binning_metabat/coassembly/{wildcards.sample_type}/*.fa"))
    conda:
        "drep"
    log:
        "logs/dereplicate_genomes/coassembly/{sample_type}.log"
    threads: 8
    resources:
        qos="normal",
        mem_mb=32000,
        time=120
    shell:
        """
        mkdir -p {params.dir_out}/genomes

        for bin in {params.bins}; do
            bin_name=$(basename "$bin" .fa)
            cp "$bin" {params.dir_out}/genomes/${{bin_name}}.fasta

        done

        [ -d {params.dir_out}/data_tables ] && rm -rf {params.dir_out}/data_tables

        dRep dereplicate {params.dir_out} \
            -g {params.dir_out}/genomes/*.fasta \
            --genomeInfo {input.checkm2_table} \
            -p {threads} \
            > {log} 2>&1

        rm -r {params.dir_out}/genomes
        """
