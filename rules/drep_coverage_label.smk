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
rule coverm_compatibility_mag:
    input:
        genomes=glob.glob("snakestream/dereplicated_genome/all/dereplicated_genomes/*.fasta"),
    output:
        dir_out=directory("snakestream/relative_abundances/dereplicated_genomes")
    conda:
        "coverm"
    benchmark:
       "benchmarks/coverm_compatibility_mag/bench.txt"
    log:
        "logs/coverm_compatibility_mag/log.log"
    threads: 8
    resources:
        qos="normal",
        mem_mb=32000,
        time=1430,
        **default_resources()
    shell:
        """
        mkdir -p {output.dir_out}
        for bin in {input.genomes}; do
            bin_name=$(basename "$bin" .fasta)
            cp "$bin" {output.dir_out}/${{bin_name}}.fna

        done
        """
#----------------------------------------------------------------------------------------#
rule coverm_relative_abundances:
    input:
        genomes_dir="snakestream/relative_abundances/dereplicated_genomes",
        r1_clean = "snakestream/reads_clean/{sample}_R1_clean.fastq.gz",
        r2_clean = "snakestream/reads_clean/{sample}_R2_clean.fastq.gz",
    output:
        abundances="snakestream/relative_abundances/{sample}_output_coverm.tsv"
    conda:
        "coverm"
    benchmark:
       "benchmarks/coverm/{sample}.txt"
    log:
        "logs/coverm_relative_abundances/{sample}.log"
    threads: 8
    resources:
        qos="normal",
        mem_mb=32000,
        time=1430,
        **default_resources()
    shell:
        """
        coverm genome \
        -1 {input.r1_clean} \
        -2 {input.r2_clean} \
        --genome-fasta-directory {input.genomes_dir} \
        -t {threads} \
        --min-covered-fraction 0\
        -m relative_abundance mean trimmed_mean covered_fraction covered_bases variance length count reads_per_base rpkm tpm\
        -o {output.abundances}
        """
#----------------------------------------------------------------------------------------#
rule gtdbtk_labeling:
    input:
        genomes_dir="snakestream/relative_abundances/dereplicated_genomes",
    output:
        log="snakestream/gtdbtk/gtdbtk.warnings.log"
    conda:
        "gtdbtk"
    params:
        out_dir="snakestream/gtdbtk",
        data_dir="PATH TO DATABASE",
    benchmark:
       "benchmarks/gtdbtk/benchmark.txt"
    log:
        "logs/gtdbtk/log.log"
    threads: 64
    resources:
        qos="normal",
        mem_mb=32000,
        time=1430,
        **default_resources()
    shell:
       """
        gtdntk classify_wf \
        --genome_dir {input.genomes_dir} \
        --out_dir {params.out_dir} \
        --data_dir {params.data_dir} \
        --cpus {threads}
       """
#----------------------------------------------------------------------------------------#
rule combine_coverm_metrics:
    input:
        expand("snakestream/relative_abundances/{sample}_output_coverm.tsv", sample=SAMPLES)
    output:
        "snakestream/coverm_combined/relative_abundance.tsv",
        "snakestream/coverm_combined/mean.tsv",
        "snakestream/coverm_combined/trimmed_mean.tsv",
        "snakestream/coverm_combined/covered_fraction.tsv",
        "snakestream/coverm_combined/covered_bases.tsv",
        "snakestream/coverm_combined/variance.tsv",
        "snakestream/coverm_combined/length.tsv",
        "snakestream/coverm_combined/read_count.tsv",
        "snakestream/coverm_combined/reads_per_base.tsv",
        "snakestream/coverm_combined/rpkm.tsv",
        "snakestream/coverm_combined/tpm.tsv"
    benchmark:
       "benchmarks/coverm_combination/benchmark.txt"
    log:
        "logs/coverm_combination/log.log"
    threads: 64
    resources:
        qos="normal",
        mem_mb=32000,
        time=1430,
        **default_resources()
    run:
        import pandas as pd, re, os
        
        metrics = {
            "Relative Abundance (%)": "relative_abundance.tsv",
            "Mean": "mean.tsv",
            "Trimmed Mean": "trimmed_mean.tsv",
            "Covered Fraction": "covered_fraction.tsv",
            "Covered Bases": "covered_bases.tsv",
            "Variance": "variance.tsv",
            "Length": "length.tsv",
            "Read Count": "read_count.tsv",
            "Reads per base": "reads_per_base.tsv",
            "RPKM": "rpkm.tsv",
            "TPM": "tpm.tsv",
        }

        dfs = {m: [] for m in metrics}

        for f in input:
            sample = re.search(r'(\d+)_R1_clean', f)
            sample = sample.group(1) if sample else os.path.basename(f).replace('.tsv', '')

            # 🔹 CoverM usa spazi multipli, li convertiamo in tab
            with open(f) as fin:
                lines = [re.sub(r'\s{2,}', '\t', line.strip()) for line in fin]
            clean_text = "\n".join(lines)

            # Leggi come TSV pulito
            from io import StringIO
            df = pd.read_csv(StringIO(clean_text), sep="\t")
            df.columns = [c.strip() for c in df.columns]

            for metric, outname in metrics.items():
                col = [c for c in df.columns if metric in c]
                if not col:
                    continue
                temp = df[['Genome', col[0]]].copy()
                temp.columns = ['Genome', sample]
                temp[sample] = pd.to_numeric(temp[sample], errors='coerce')
                dfs[metric].append(temp)

        os.makedirs("snakestream/coverm_combined", exist_ok=True)

        for metric, outname in metrics.items():
            if not dfs[metric]:
                continue
            merged = dfs[metric][0]
            for d in dfs[metric][1:]:
                merged = merged.merge(d, on='Genome', how='outer')
            merged.to_csv(f"snakestream/coverm_combined/{outname}", sep='\t', index=False)
            pd.read_csv(f"snakestream/coverm_combined/{outname}",sep=r'\s{2,}|\t', engine='python').to_csv(f"snakestream/coverm_combined/{outname}",sep="\t",index=False)
        sys.stdout = open(log[0], "w")
        sys.stderr = sys.stdout
