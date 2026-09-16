# Raw read statistics.
rule seqkit_raw_reads:
    input: expand("data/reads_raw/{sample}_{read}.fastq.gz", sample=SAMPLES, read=["R1", "R2"])
    output: protected("tables/read_stats/seqkit_raw_reads.tsv")
    log: "logs/stats/seqkit_raw_reads.log"
    benchmark: "benchmarks/stats/seqkit_raw_reads.tsv"
    threads: 16
    resources: mem_mb=20000, runtime=600
    container: config["containers"]["seqkit"]
    shell: "seqkit stats -j {threads} -T {input:q} -o {output:q} > {log:q} 2>&1"

# Trimmed read statistics.
rule seqkit_trimmed_reads:
    input: expand("data/reads_trim/{sample}_{read}_trim.fastq.gz", sample=SAMPLES, read=["R1", "R2"])
    output: protected("tables/read_stats/seqkit_trimmed_reads.tsv")
    log: "logs/stats/seqkit_trimmed_reads.log"
    benchmark: "benchmarks/stats/seqkit_trimmed_reads.tsv"
    threads: 16
    resources: mem_mb=20000, runtime=600
    container: config["containers"]["seqkit"]
    shell: "seqkit stats -j {threads} -T {input:q} -o {output:q} > {log:q} 2>&1"

# Cleaned read statistics.
rule seqkit_cleaned_reads:
    input: expand("data/reads_clean/{sample}_{read}_clean.fastq.gz", sample=SAMPLES, read=["R1", "R2"])
    output: protected("tables/read_stats/seqkit_cleaned_reads.tsv")
    log: "logs/stats/seqkit_cleaned_reads.log"
    benchmark: "benchmarks/stats/seqkit_cleaned_reads.tsv"
    threads: 16
    resources: mem_mb=20000, runtime=600
    container: config["containers"]["seqkit"]
    shell: "seqkit stats -j {threads} -T {input:q} -o {output:q} > {log:q} 2>&1"

# Singleton read statistics.
rule seqkit_cleaned_reads_sing:
    input: expand("data/reads_clean/{sample}_sing_clean.fastq.gz", sample=SAMPLES)
    output: protected("tables/read_stats/seqkit_cleaned_reads_sing.tsv")
    log: "logs/stats/seqkit_cleaned_reads_sing.log"
    benchmark: "benchmarks/stats/seqkit_cleaned_reads_sing.tsv"
    threads: 16
    resources: mem_mb=20000, runtime=600
    container: config["containers"]["seqkit"]
    shell: "seqkit stats -j {threads} -T {input:q} -o {output:q} > {log:q} 2>&1"
