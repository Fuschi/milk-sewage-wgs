HOST_INDEX = config["references"]["bostaurus_index"]
HOST_INDEX_FILES = expand(HOST_INDEX + ".{ext}", ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])

# Trim reads.
rule bbduk_trim:
    input:
        r1="data/reads_raw/{sample}_R1.fastq.gz",
        r2="data/reads_raw/{sample}_R2.fastq.gz",
    output:
        r1=protected("data/reads_trim/{sample}_R1_trim.fastq.gz"),
        r2=protected("data/reads_trim/{sample}_R2_trim.fastq.gz"),
        singleton=protected("data/reads_trim/{sample}_sing.fastq.gz"),
        stats=protected("data/reads_trim/{sample}_stats.txt"),
    params: adapters=config["references"]["adapters"]
    log: "logs/bbduk_trim/{sample}.log"
    benchmark: "benchmarks/bbduk_trim/{sample}.tsv"
    threads: 4
    resources: mem_mb=15000, runtime=30
    container: config["containers"]["bbmap"]
    shell:
        """
        bbduk.sh -Xmx14g threads={threads} \
            in1={input.r1} in2={input.r2} \
            out1={output.r1} out2={output.r2} \
            outs={output.singleton} stats={output.stats} \
            ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tpe tbo \
            qtrim=rl trimq=10 ow=t ziplevel=6 > {log} 2>&1
        """

# Filter host reads from paired-end data.
rule remove_host:
    input:
        r1="data/reads_trim/{sample}_R1_trim.fastq.gz",
        r2="data/reads_trim/{sample}_R2_trim.fastq.gz",
        index=HOST_INDEX_FILES,
    output:
        r1=protected("data/reads_clean/{sample}_R1_clean.fastq.gz"),
        r2=protected("data/reads_clean/{sample}_R2_clean.fastq.gz"),
        bam=temp("data/reads_clean/{sample}_hostaligned.bam"),
        unmapped=temp("data/reads_clean/{sample}_unmapped.bam"),
        sorted=protected("data/reads_clean/{sample}_unmapped_sorted.bam"),
    params: index=HOST_INDEX
    log: "logs/remove_host/{sample}.log"
    benchmark: "benchmarks/remove_host/{sample}.tsv"
    threads: 4
    resources: mem_mb=64000, runtime=180
    container: config["containers"]["hostremoval"]
    shell:
        """
        (
        bowtie2 -p {threads} -x {params.index} -1 {input.r1} -2 {input.r2} -S /dev/stdout \
            | samtools view -bS - > {output.bam}
        samtools view -b -f 12 -F 256 {output.bam} > {output.unmapped}
        samtools sort -n -m 5G -@ {threads} {output.unmapped} -o {output.sorted}
        samtools fastq -@ {threads} -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {output.sorted}
        ) > {log} 2>&1
        """

# Filter host reads from single-end data.
rule remove_host_sing:
    input:
        singleton="data/reads_trim/{sample}_sing.fastq.gz",
        index=HOST_INDEX_FILES,
    output:
        singleton=protected("data/reads_clean/{sample}_sing_clean.fastq.gz"),
        bam=temp("data/reads_clean/{sample}_sing_hostaligned.bam"),
        unmapped=temp("data/reads_clean/{sample}_sing_unmapped.bam"),
        sorted=protected("data/reads_clean/{sample}_sing_unmapped_sorted.bam"),
    params: index=HOST_INDEX
    log: "logs/remove_host/{sample}_sing.log"
    benchmark: "benchmarks/remove_host/{sample}_sing.tsv"
    threads: 4
    resources: mem_mb=64000, runtime=180
    container: config["containers"]["hostremoval"]
    shell:
        """
        (
        bowtie2 -p {threads} -x {params.index} -U {input.singleton} -S /dev/stdout \
            | samtools view -bS - > {output.bam}
        samtools view -b -f 12 -F 256 {output.bam} > {output.unmapped}
        samtools sort -n -m 5G -@ {threads} {output.unmapped} -o {output.sorted}
        samtools fastq -@ {threads} -0 {output.singleton} -s /dev/null -n {output.sorted}
        ) > {log} 2>&1
        """
