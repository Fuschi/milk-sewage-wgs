MATRIX_PARTITIONS = "m1,m2,m3,m4,m5,m6,m7"
COVERM_METHODS = "mean trimmed_mean covered_fraction covered_bases variance length count reads_per_base anir rpkm tpm"

# Concatenate individual assemblies by biome for multi-sample binning.
rule multisample_fasta:
    input: lambda wc: expand("data/assembly/sample/{sample}/{sample}.contigs.fa", sample=SAMPLES_BY_BIOME[wc.biome])
    output: "data/assembly/multisample/{biome}/{biome}.contigs.fa"
    resources: mem_mb=4000, runtime=30
    shell:
        """
        cat {input} > {output}
        """

# Build Strobealign index for each individual assembly.
rule strobealign_index_sample:
    input: ref="data/assembly/sample/{sample}/{sample}.contigs.fa"
    output: "data/assembly/sample/{sample}/{sample}.contigs.fa.r150.sti"
    threads: 8
    resources: mem_mb=32000, runtime=60, slurm_partition=MATRIX_PARTITIONS
    container: config["containers"]["strobealign"]
    log: "logs/contig_mapping/sample/{sample}_index.log"
    benchmark: "benchmarks/contig_mapping/sample/{sample}_index.tsv"
    shell:
        """
        strobealign --create-index -r 150 -t {threads} {input.ref} > {log} 2>&1
        """


# Map paired and singleton reads to the individual assembly.
rule contig_mapping_sample:
    input:
        ref="data/assembly/sample/{sample}/{sample}.contigs.fa",
        index="data/assembly/sample/{sample}/{sample}.contigs.fa.r150.sti",
        r1="data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz",
        r2="data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz",
        sing1="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz",
        sing2="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz"
    output:
        bam="data/contig_mapping/sample/{sample}/{sample}.sorted.bam",
        bai="data/contig_mapping/sample/{sample}/{sample}.sorted.bam.bai"
    threads: 16
    resources: mem_mb=32000, runtime=240, slurm_partition=MATRIX_PARTITIONS
    container: config["containers"]["strobealign"]
    log: "logs/contig_mapping/sample/{sample}.log"
    benchmark: "benchmarks/contig_mapping/sample/{sample}.tsv"
    shell:
        """
        strobealign --use-index -r 150 -U -t 12 {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort -@ 4 -m 3G -o {output.bam}.pe.bam -
        strobealign --use-index -r 150 -U -t 12 {input.ref} {input.sing1} 2>> {log} | samtools sort -@ 4 -m 3G -o {output.bam}.s1.bam -
        strobealign --use-index -r 150 -U -t 12 {input.ref} {input.sing2} 2>> {log} | samtools sort -@ 4 -m 3G -o {output.bam}.s2.bam -
        samtools merge -@ 4 -f {output.bam} {output.bam}.pe.bam {output.bam}.s1.bam {output.bam}.s2.bam 2>> {log}
        samtools index -@ 4 {output.bam} {output.bai} 2>> {log}
        rm {output.bam}.pe.bam {output.bam}.s1.bam {output.bam}.s2.bam
        """


# Calculate coverage metrics from the individual assembly BAM.
rule contig_coverage_sample:
    input: "data/contig_mapping/sample/{sample}/{sample}.sorted.bam"
    output: "data/contig_mapping/sample/{sample}/{sample}.coverage.tsv"
    params: methods=COVERM_METHODS
    threads: 4
    resources: mem_mb=8000, runtime=30
    container: config["containers"]["coverm"]
    log: "logs/contig_mapping/sample/{sample}_coverage.log"
    benchmark: "benchmarks/contig_mapping/sample/{sample}_coverage.tsv"
    shell:
        """
        coverm contig --bam-files {input} --methods {params.methods} -t {threads} -o {output} > {log} 2>&1
        """


# Build Strobealign index for each biome coassembly.
rule strobealign_index_coassembly:
    input: ref="data/assembly/coassembly/{biome}/{biome}.contigs.fa"
    output: "data/assembly/coassembly/{biome}/{biome}.contigs.fa.r150.sti"
    threads: 16
    resources: mem_mb=128000, runtime=120, slurm_partition=MATRIX_PARTITIONS
    container: config["containers"]["strobealign"]
    log: "logs/contig_mapping/coassembly/{biome}_index.log"
    benchmark: "benchmarks/contig_mapping/coassembly/{biome}_index.tsv"
    shell:
        """
        strobealign --create-index -r 150 -t {threads} {input.ref} > {log} 2>&1
        """


# Map paired and singleton reads to the biome coassembly.
rule contig_mapping_coassembly:
    input:
        ref="data/assembly/coassembly/{biome}/{biome}.contigs.fa",
        index="data/assembly/coassembly/{biome}/{biome}.contigs.fa.r150.sti",
        r1="data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz",
        r2="data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz",
        sing1="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz",
        sing2="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz"
    output:
        bam="data/contig_mapping/coassembly/{biome}/{sample}.sorted.bam",
        bai="data/contig_mapping/coassembly/{biome}/{sample}.sorted.bam.bai"
    threads: 16
    resources: mem_mb=128000, runtime=360, slurm_partition=MATRIX_PARTITIONS
    container: config["containers"]["strobealign"]
    log: "logs/contig_mapping/coassembly/{biome}_{sample}.log"
    benchmark: "benchmarks/contig_mapping/coassembly/{biome}_{sample}.tsv"
    shell:
        """
        strobealign --use-index -r 150 -U -t 12 {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort -@ 4 -m 3G -o {output.bam}.pe.bam -
        strobealign --use-index -r 150 -U -t 12 {input.ref} {input.sing1} 2>> {log} | samtools sort -@ 4 -m 3G -o {output.bam}.s1.bam -
        strobealign --use-index -r 150 -U -t 12 {input.ref} {input.sing2} 2>> {log} | samtools sort -@ 4 -m 3G -o {output.bam}.s2.bam -
        samtools merge -@ 4 -f {output.bam} {output.bam}.pe.bam {output.bam}.s1.bam {output.bam}.s2.bam 2>> {log}
        samtools index -@ 4 {output.bam} {output.bai} 2>> {log}
        rm {output.bam}.pe.bam {output.bam}.s1.bam {output.bam}.s2.bam
        """


# Calculate coverage metrics from each coassembly BAM.
rule contig_coverage_coassembly:
    input: "data/contig_mapping/coassembly/{biome}/{sample}.sorted.bam"
    output: "data/contig_mapping/coassembly/{biome}/{sample}.coverage.tsv"
    params: methods=COVERM_METHODS
    threads: 4
    resources: mem_mb=64000, runtime=120, slurm_partition=MATRIX_PARTITIONS
    container: config["containers"]["coverm"]
    log: "logs/contig_mapping/coassembly/{biome}_{sample}_coverage.log"
    benchmark: "benchmarks/contig_mapping/coassembly/{biome}_{sample}_coverage.tsv"
    shell:
        """
        coverm contig --bam-files {input} --methods {params.methods} -t {threads} -o {output} > {log} 2>&1
        """


# Build Strobealign index for concatenated individual assemblies.
rule strobealign_index_multisample:
    input: ref="data/assembly/multisample/{biome}/{biome}.contigs.fa"
    output: "data/assembly/multisample/{biome}/{biome}.contigs.fa.r150.sti"
    threads: 16
    resources: mem_mb=128000, runtime=180, slurm_partition=MATRIX_PARTITIONS
    container: config["containers"]["strobealign"]
    log: "logs/contig_mapping/multisample/{biome}_index.log"
    benchmark: "benchmarks/contig_mapping/multisample/{biome}_index.tsv"
    shell:
        """
        strobealign --create-index -r 150 -t {threads} {input.ref} > {log} 2>&1
        """


# Map paired and singleton reads to the concatenated assemblies of the biome.
rule contig_mapping_multisample:
    input:
        ref="data/assembly/multisample/{biome}/{biome}.contigs.fa",
        index="data/assembly/multisample/{biome}/{biome}.contigs.fa.r150.sti",
        r1="data/kneaddata/{sample}/{sample}_kneaddata_paired_1.fastq.gz",
        r2="data/kneaddata/{sample}/{sample}_kneaddata_paired_2.fastq.gz",
        sing1="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_1.fastq.gz",
        sing2="data/kneaddata/{sample}/{sample}_kneaddata_unmatched_2.fastq.gz"
    output:
        bam="data/contig_mapping/multisample/{biome}/{sample}.sorted.bam",
        bai="data/contig_mapping/multisample/{biome}/{sample}.sorted.bam.bai"
    threads: 16
    resources: mem_mb=128000, runtime=480, slurm_partition=MATRIX_PARTITIONS
    container: config["containers"]["strobealign"]
    log: "logs/contig_mapping/multisample/{biome}_{sample}.log"
    benchmark: "benchmarks/contig_mapping/multisample/{biome}_{sample}.tsv"
    shell:
        """
        strobealign --use-index -r 150 -U -t 12 {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort -@ 4 -m 3G -o {output.bam}.pe.bam -
        strobealign --use-index -r 150 -U -t 12 {input.ref} {input.sing1} 2>> {log} | samtools sort -@ 4 -m 3G -o {output.bam}.s1.bam -
        strobealign --use-index -r 150 -U -t 12 {input.ref} {input.sing2} 2>> {log} | samtools sort -@ 4 -m 3G -o {output.bam}.s2.bam -
        samtools merge -@ 4 -f {output.bam} {output.bam}.pe.bam {output.bam}.s1.bam {output.bam}.s2.bam 2>> {log}
        samtools index -@ 4 {output.bam} {output.bai} 2>> {log}
        rm {output.bam}.pe.bam {output.bam}.s1.bam {output.bam}.s2.bam
        """


# Calculate coverage metrics from each multi-sample BAM.
rule contig_coverage_multisample:
    input: "data/contig_mapping/multisample/{biome}/{sample}.sorted.bam"
    output: "data/contig_mapping/multisample/{biome}/{sample}.coverage.tsv"
    params: methods=COVERM_METHODS
    threads: 4
    resources: mem_mb=128000, runtime=180, slurm_partition=MATRIX_PARTITIONS
    container: config["containers"]["coverm"]
    log: "logs/contig_mapping/multisample/{biome}_{sample}_coverage.log"
    benchmark: "benchmarks/contig_mapping/multisample/{biome}_{sample}_coverage.tsv"
    shell:
        """
        coverm contig --bam-files {input} --methods {params.methods} -t {threads} -o {output} > {log} 2>&1
        """

