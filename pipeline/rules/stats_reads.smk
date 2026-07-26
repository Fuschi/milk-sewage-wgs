# -----------------------------------------------------------------------------
# Read QC and summary statistics
# -----------------------------------------------------------------------------

rule fastqc_raw:
    input:
        r1="reads_raw/{sample}_R1.fastq.gz",
        r2="reads_raw/{sample}_R2.fastq.gz",
    output:
        html_r1=protected("qc/raw/{sample}_R1_fastqc.html"),
        html_r2=protected("qc/raw/{sample}_R2_fastqc.html"),
    log:
        "logs/qc/raw/{sample}.log"
    threads: 2
    resources:
        mem_mb=5000,
        time="00:15:00",
        qos="normal",
    conda:
        str(ENVS_DIR / "fastqc.yml")
    shell:
        r"""
        set -euo pipefail
        mkdir -p qc/raw logs/qc/raw
        fastqc {input.r1} {input.r2} \
            --threads {threads} \
            --outdir qc/raw \
            > {log} 2>&1
        """


rule fastqc_trimmed:
    input:
        r1="reads_trim/{sample}_R1_trim.fastq.gz",
        r2="reads_trim/{sample}_R2_trim.fastq.gz",
    output:
        html_r1=protected("qc/trim/{sample}_R1_trim_fastqc.html"),
        html_r2=protected("qc/trim/{sample}_R2_trim_fastqc.html"),
    log:
        "logs/qc/trim/{sample}.log"
    threads: 2
    resources:
        mem_mb=2000,
        time="00:10:00",
        qos="normal",
    conda:
        str(ENVS_DIR / "fastqc.yml")
    shell:
        r"""
        set -euo pipefail
        mkdir -p qc/trim logs/qc/trim
        fastqc {input.r1} {input.r2} \
            --threads {threads} \
            --outdir qc/trim \
            > {log} 2>&1
        """


rule seqkit_raw_reads:
    input:
        expand(
            "reads_raw/{sample}_{pair}.fastq.gz",
            sample=SAMPLES,
            pair=["R1", "R2"],
        )
    output:
        protected(str(READ_STATS_DIR / "seqkit_raw_reads.tsv"))
    threads: 16
    resources:
        mem_mb=20000,
        time="10:00:00",
        qos="normal",
    conda:
        str(ENVS_DIR / "seqkit.yml")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        seqkit stats -j {threads} -T {input} -o {output}
        """


rule seqkit_trimmed_reads:
    input:
        expand(
            "reads_trim/{sample}_{pair}_trim.fastq.gz",
            sample=SAMPLES,
            pair=["R1", "R2"],
        )
    output:
        protected(str(READ_STATS_DIR / "seqkit_trimmed_reads.tsv"))
    threads: 16
    resources:
        mem_mb=20000,
        time="10:00:00",
        qos="normal",
    conda:
        str(ENVS_DIR / "seqkit.yml")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        seqkit stats -j {threads} -T {input} -o {output}
        """


rule seqkit_cleaned_reads:
    input:
        expand(
            "reads_clean/{sample}_{pair}_clean.fastq.gz",
            sample=SAMPLES,
            pair=["R1", "R2"],
        )
    output:
        protected(str(READ_STATS_DIR / "seqkit_cleaned_reads.tsv"))
    threads: 16
    resources:
        mem_mb=20000,
        time="10:00:00",
        qos="normal",
    conda:
        str(ENVS_DIR / "seqkit.yml")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        seqkit stats -j {threads} -T {input} -o {output}
        """
