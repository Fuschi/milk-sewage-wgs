# rules/cleaning.smk
# ─────────────────────────────────────────────
# Filtering and cleaning rules:
# - bostaurus_reference
# - bbduk_trim
# - remove_host
# ─────────────────────────────────────────────
#Generic Resources
def default_resources():
    return dict(
        requeue=0,
        trigger=1,
    )
# Download Bos taurus host genome from UCSC and build Bowtie2 index
# ─────────────────────────────────────────────────────────────
rule bostaurus_reference:
    output:
        fastagz = protected("snakestream/host/bosTau9.fa.gz"),
        fasta   = temp("snakestream/host/bosTau9.fa"),
        index_files = protected(expand("snakestream/host/bosTau9.{ext}", ext=[
            "1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"
        ]))
    log:
        "logs/bostaurus_reference.log"
    threads: 4
    resources:
        mem_mb=10000,
        time="01:00:00",
        qos="normal",
        **default_resources(),
    conda:
        "hostremoval"
    shell:
        """
        (
        set -euo pipefail

        echo "== Downloading bosTau9 reference genome =="
        wget -O {output.fastagz} \
            https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/bigZips/bosTau9.fa.gz 

        echo ""; echo "== Decompressing FASTA =="
        gunzip -c {output.fastagz} > {output.fasta} 

        echo ""; echo "== Building Bowtie2 index =="
        bowtie2-build --threads {threads} --seed 42 {output.fasta} snakestream/host/bosTau9
        echo "== Done =="
		) > {log} 2>&1
        """

# Trim raw reads and remove adapters using BBDuk
# ─────────────────────────────────────────────────────────────
rule bbduk_trim:
    input:
        r1 = "snakestream/reads_raw/{sample}_R1.fastq.gz",
        r2 = "snakestream/reads_raw/{sample}_R2.fastq.gz",
        adapters = "/home/PERSONALE/alessandro.fuschi2/miniforge3/pkgs/bbmap-36.32-0/opt/bbmap-36.32/resources/adapters.fa"
    output:
        r1_trim = protected("snakestream/reads_trim/{sample}_R1_trim.fastq.gz"),
        r2_trim = protected("snakestream/reads_trim/{sample}_R2_trim.fastq.gz"),
        singleton = protected("snakestream/reads_trim/{sample}_sing.fastq.gz"),
        stats = protected("snakestream/reads_trim/{sample}_stats.txt")
    log:
        "logs/bbduk_trim/{sample}.log"
    threads: 4
    resources:
        mem_mb=15000,
        time="00:30:00",
        qos="normal",
        **default_resources(),
    conda:
        "bbmap"
    shell:
        """
		bbduk.sh \
            -Xmx14g \
            in1={input.r1} in2={input.r2} \
            out1={output.r1_trim} out2={output.r2_trim} \
            outs={output.singleton}\
            stats={output.stats} \
            ref={input.adapters} \
            ktrim=r k=23 mink=11 hdist=1 \
            tpe tbo \
            qtrim=rl trimq=10 \
            ow=t ziplevel=6 \
			> {log} 2>&1
        """

# Remove host reads from R1 and R2 using Bowtie2 and Samtools
# ─────────────────────────────────────────────────────────────
rule remove_host:
    input:
        r1 = "snakestream/reads_trim/{sample}_R1_trim.fastq.gz",
        r2 = "snakestream/reads_trim/{sample}_R2_trim.fastq.gz",
        index_check = "snakestream/host/bosTau9.1.bt2"
    output:
        r1_clean = protected("snakestream/reads_clean/{sample}_R1_clean.fastq.gz"),
        r2_clean = protected("snakestream/reads_clean/{sample}_R2_clean.fastq.gz"),
        bam =      temp("snakestream/reads_clean/{sample}_hostaligned.bam"),
        unmapped = temp("snakestream/reads_clean/{sample}_unmapped.bam"),
        sorted =   protected("snakestream/reads_clean/{sample}_unmapped_sorted.bam")
    params:
        index = "snakestream/host/bosTau9"
    log:
        "logs/remove_host/{sample}.log"
    threads: 4
    resources:
        mem_mb=64000,
        time="03:00:00",
        qos="normal",
        **default_resources(),
    conda:
        "hostremoval"
    shell:
        """
        (
        set -euo pipefail

        echo "== Bowtie2 mapping =="
        bowtie2 -p {threads} -x {params.index} \
             -1 {input.r1} -2 {input.r2} -S /dev/stdout |
        samtools view -bS - > {output.bam}

        echo "== Filtering unmapped reads =="
        samtools view -b -f 12 -F 256 {output.bam} > {output.unmapped}

        echo "== Sorting BAM =="
        samtools sort -n -m 5G -@ {threads} {output.unmapped} -o {output.sorted}

        echo "== Converting BAM to FASTQ =="
        samtools fastq -@ {threads} -1 {output.r1_clean} -2 {output.r2_clean} \
        -0 /dev/null -s /dev/null -n {output.sorted}

        echo "== Done =="
        ) > {log} 2>&1
        """
# Remove host reads from singleton using Bowtie2 and Samtools
# ─────────────────────────────────────────────────────────────
rule remove_host_sing:
    input:
        singleton = "snakestream/reads_trim/{sample}_sing.fastq.gz",
        index_check = "snakestream/host/bosTau9.1.bt2"
    output:
        singleton_clean = protected("snakestream/reads_clean/{sample}_sing_clean.fastq.gz"),
        bam =      temp("snakestream/reads_clean/{sample}_sing_hostaligned.bam"),
        unmapped = temp("snakestream/reads_clean/{sample}_sing_unmapped.bam"),
        sorted =   protected("snakestream/reads_clean/{sample}_sing_unmapped_sorted.bam")
    params:
        index = "snakestream/host/bosTau9"
    log:
        "logs/remove_host/{sample}_sing.log"
    threads: 4
    resources:
        mem_mb=64000,
        time="03:00:00",
        qos="normal",
        **default_resources(),
    conda:
        "hostremoval"
    shell:
        """
        (
        set -euo pipefail

        echo "== Bowtie2 mapping =="
        bowtie2 -p {threads} -x {params.index} \
             -U {input.singleton} -S /dev/stdout |
        samtools view -bS - > {output.bam}

        echo "== Filtering unmapped reads =="
        samtools view -b -f 12 -F 256 {output.bam} > {output.unmapped}

        echo "== Sorting BAM =="
        samtools sort -n -m 5G -@ {threads} {output.unmapped} -o {output.sorted}

        echo "== Converting BAM to FASTQ =="
        samtools fastq -@ {threads} \
         -0 {output.singleton_clean} -s /dev/null -n {output.sorted}

        echo "== Done =="
        ) > {log} 2>&1
        """
