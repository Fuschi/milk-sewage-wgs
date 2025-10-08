# rules/drep_coverage_label.smk
# ─────────────────────────────────────────────
import glob
# Binning for assembly conting rules:
# ─────────────────────────────────────────────
#Generic Resources
def default_resources():
    return dict(
        requeue=0,
        trigger=1,
    )
#----------------------------------------------------------------------------------------#
rule merge_checkm2_reports_all:
    input:
        log = [glob.glob("snakestream/checkm2/assembly/*/*/checkm2.log")+glob.glob("snakestream/checkm2/coassembly/*/checkm2.log")]
    output:
        tables="snakestream/tables/all/checkm2_reports_for_dRep.csv"
    params:
        reports=[glob.glob("snakestream/checkm2/assembly/*/*/*quality_report.tsv")+glob.glob("snakestream/checkm2/coassembly/*/*quality_report.tsv")] 
    benchmark:
       "benchmarks/merge_checkm2/all/all.txt"
    resources:
        qos="normal",
        mem_mb=32000,
        time=120,
        **default_resources()
    conda:
        "r-tidyverse"
    script:
        "../scripts/merge_checkm2_reports_assembly.R"


#----------------------------------------------------------------------------------------#
rule dereplicate_genomes_all:
    input:
          bins=[glob.glob("snakestream/binning_metabat/*/*/*_bin.BinInfo.txt")+glob.glob("snakestream/binning_metabat/*/*/*/*_bin.BinInfo.txt")],
          checkm2_table="snakestream/tables/all/checkm2_reports_for_dRep.csv"
    output:
        "snakestream/dereplicated_genome/all/figures/Winning_genomes.pdf"
    params:
        dir_out="snakestream/dereplicated_genome/all",
        bins=glob.glob("snakestream/binning_metabat/*/*/*.fa")+glob.glob("snakestream/binning_metabat/*/*/*/*.fa")
    benchmark:
       "benchmarks/drep/all.txt"
    conda:
        "drep"
    log:
        "logs/dereplicate_genomes/all.log"
    threads: 8
    resources:
        qos="normal",
        mem_mb=32000,
        time=1430,
        **default_resources()
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
#----------------------------------------------------------------------------------------#
rule coverm_relative_abundances_assembly:
    input:
          genomes=glob.glob("snakestream/dereplicated_genome/all/dereplicated_genome/*.fasta")
          contigs=
    output:
        abundances="snakestream/relative_abundances/assembly/output_coverm_assembly.tsv"
    params:
    conda:
        "coverm"
    benchmark:
       "benchmarks/coverm/all/all.txt"
    log:
        "logs/coverm_relative_abundances/assembly/all.log"
    threads: 8
    resources:
        qos="normal",
        mem_mb=32000,
        time=1430,
        **default_resources()
    shell:
        """
        coverm genome \
        --coupled sample_1.1.fq.gz sample_1.2.fq.gz \
        --genome-fasta-files {input.genomes} \
        -t {threads} \
        -m relative_abundance  \
        -o {output.abundances}
        """
