# rules/reads_stats.smk
# ─────────────────────────────────────────────
# Getting stats from reads:
# - fastqc_raw / fastqc_trimmed
# - seqkit_raw_reads / trim_reads / clean_reads
# ─────────────────────────────────────────────
#Generic Resources
def default_resources():
    return dict(
        requeue=0,
        trigger=1,
    )
# Raw reads QC and stats
rule fastqc_raw:
    input:
        r1 = "snakestream/reads_raw/{sample}_R1.fastq.gz",
        r2 = "snakestream/reads_raw/{sample}_R2.fastq.gz"
    output:
        html_r1 = protected("snakestream/qc/raw/{sample}_R1_fastqc.html"),
        html_r2 = protected("snakestream/qc/raw/{sample}_R2_fastqc.html")
    log:
        "logs/qc/raw/{sample}.log"
    threads: 2
    resources:
        mem_mb=5000, 
        time="00:15:00", 
        qos="normal",
        **default_resources(),
    conda:
        "fastqc"
    shell:
        "fastqc {input.r1} {input.r2} -t {threads} -o snakestream/qc/raw > {log} 2>&1"

rule seqkit_raw_reads:
    input:
        expand("snakestream/reads_raw/{sample}_{pe}.fastq.gz", sample=SAMPLES, pe=["R1", "R2"])
    output:
        protected("snakestream/stats/seqkit_raw_reads.tsv")
    threads: 16
    resources:
        mem_mb=20000, 
        time="10:00:00", 
        qos="normal",
        **default_resources(),
    conda:
        "seqkit"
    shell:
        "seqkit stats -j {threads} -T {input} -o {output}"

# Trimmed reads QC and stats
rule fastqc_trimmed:
    input:
        r1 = "snakestream/reads_trim/{sample}_R1_trim.fastq.gz",
        r2 = "snakestream/reads_trim/{sample}_R2_trim.fastq.gz"
    output:
        html_r1 = protected("snakestream/qc/trim/{sample}_R1_trim_fastqc.html"),
        html_r2 = protected("snakestream/qc/trim/{sample}_R2_trim_fastqc.html")
    log:
        "logs/qc/trim/{sample}.log"
    threads: 2
    resources:
        mem_mb=2000, 
        time="00:10:00", 
        qos="normal",
        **default_resources(),
    conda:
        "fastqc"
    shell:
        "fastqc {input.r1} {input.r2} -t {threads} -o snakestream/qc/trim > {log} 2>&1"

rule seqkit_trimmed_reads:
    input:
        expand("snakestream/reads_trim/{sample}_{pe}_trim.fastq.gz", sample=SAMPLES, pe=["R1", "R2"])
    output:
        protected("snakestream/stats/seqkit_trimmed_reads.tsv")
    threads: 16
    resources:
        mem_mb=20000,
        time="10:00:00",
        qos="normal",
        **default_resources()
    conda:
        "seqkit"
    shell:
        "seqkit stats -j {threads} -T {input} -o {output}"

# Cleaned reads stats
rule seqkit_cleaned_reads:
    input:
        expand("snakestream/reads_clean/{sample}_{pe}_clean.fastq.gz", sample=SAMPLES, pe=["R1", "R2"])
    output:
        protected("snakestream/stats/seqkit_cleaned_reads.tsv")
    threads: 16
    resources:
        mem_mb=20000, 
        time="10:00:00", 
        qos="normal",
        **default_resources(),
    conda:
        "seqkit"
    shell:
        "seqkit stats -j {threads} -T {input} -o {output}"
# Cleaned sing reads stats
rule seqkit_cleaned_reads_sing:
    input:
        expand("snakestream/reads_clean/{sample}_sing_clean.fastq.gz", sample=SAMPLES)
    output:
        protected("snakestream/stats/seqkit_cleaned_reads_sing.tsv")
    threads: 16
    resources:
        mem_mb=20000,
        time="10:00:00",
        qos="normal",
        **default_resources(),
    conda:
        "seqkit"
    shell:
        "seqkit stats -j {threads} -T {input} -o {output}"
