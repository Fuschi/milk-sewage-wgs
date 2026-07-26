# -----------------------------------------------------------------------------
# Read cleaning
# -----------------------------------------------------------------------------

rule bostaurus_reference:
    output:
        fastagz=protected("host/bosTau9.fa.gz"),
        fasta=temp("host/bosTau9.fa"),
        index_files=protected(
            expand(
                "host/bosTau9.{ext}",
                ext=[
                    "1.bt2",
                    "2.bt2",
                    "3.bt2",
                    "4.bt2",
                    "rev.1.bt2",
                    "rev.2.bt2",
                ],
            )
        ),
    log:
        "logs/bostaurus_reference.log"
    threads: 4
    resources:
        mem_mb=10000,
        time="01:00:00",
        qos="normal",
    conda:
        str(ENVS_DIR / "hostremoval.yml")
    shell:
        r"""
        set -euo pipefail
        mkdir -p host logs

        (
        echo "== Downloading bosTau9 reference genome =="
        wget -q -O {output.fastagz} \
            https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/bigZips/bosTau9.fa.gz

        echo "== Decompressing FASTA =="
        gunzip -c {output.fastagz} > {output.fasta}

        echo "== Building Bowtie2 index =="
        bowtie2-build --threads {threads} --seed 42 \
            {output.fasta} host/bosTau9

        echo "== Done =="
        ) > {log} 2>&1
        """


rule bbduk_trim:
    input:
        r1="reads_raw/{sample}_R1.fastq.gz",
        r2="reads_raw/{sample}_R2.fastq.gz",
    output:
        r1_trim=protected("reads_trim/{sample}_R1_trim.fastq.gz"),
        r2_trim=protected("reads_trim/{sample}_R2_trim.fastq.gz"),
        singleton=protected("reads_trim/{sample}_sing.fastq.gz"),
        stats=protected("reads_trim/{sample}_stats.txt"),
    params:
        # Default adapter file bundled with the pinned BBMap package.
        # Override with: --config bbmap_adapters=/path/adapters.fa
        adapters=config.get(
            "bbmap_adapters",
            "$CONDA_PREFIX/opt/bbmap-36.32/resources/adapters.fa",
        ),
    log:
        "logs/bbduk_trim/{sample}.log"
    threads: 4
    resources:
        mem_mb=15000,
        time="00:30:00",
        qos="normal",
    conda:
        str(ENVS_DIR / "bbmap.yml")
    shell:
        r"""
        set -euo pipefail
        mkdir -p reads_trim logs/bbduk_trim

        bbduk.sh \
            in1={input.r1} \
            in2={input.r2} \
            out1={output.r1_trim} \
            out2={output.r2_trim} \
            outs={output.singleton} \
            stats={output.stats} \
            ref="{params.adapters}" \
            ktrim=r \
            k=23 \
            mink=11 \
            hdist=1 \
            tpe \
            tbo \
            qtrim=rl \
            trimq=10 \
            ow=t \
            ziplevel=6 \
            > {log} 2>&1
        """


rule remove_host:
    input:
        r1="reads_trim/{sample}_R1_trim.fastq.gz",
        r2="reads_trim/{sample}_R2_trim.fastq.gz",
        index_check="host/bosTau9.1.bt2",
    output:
        r1_clean=protected("reads_clean/{sample}_R1_clean.fastq.gz"),
        r2_clean=protected("reads_clean/{sample}_R2_clean.fastq.gz"),
        bam=temp("reads_clean/{sample}_hostaligned.bam"),
        unmapped=temp("reads_clean/{sample}_unmapped.bam"),
        sorted=temp("reads_clean/{sample}_unmapped_sorted.bam"),
    params:
        index="host/bosTau9",
    log:
        "logs/remove_host/{sample}.log"
    threads: 4
    resources:
        mem_mb=20000,
        time="01:00:00",
        qos="normal",
    conda:
        str(ENVS_DIR / "hostremoval.yml")
    shell:
        r"""
        set -euo pipefail
        mkdir -p reads_clean logs/remove_host

        (
        echo "== Bowtie2 mapping =="
        bowtie2 -p {threads} -x {params.index} \
            -1 {input.r1} -2 {input.r2} -S /dev/stdout \
          | samtools view -bS - \
          > {output.bam}

        echo "== Filtering unmapped read pairs =="
        samtools view -b -f 12 -F 256 \
            {output.bam} \
          > {output.unmapped}

        echo "== Sorting BAM by read name =="
        samtools sort -n -m 5G -@ {threads} \
            {output.unmapped} \
            -o {output.sorted}

        echo "== Converting BAM to paired FASTQ =="
        samtools fastq -@ {threads} \
            -1 {output.r1_clean} \
            -2 {output.r2_clean} \
            -0 /dev/null \
            -s /dev/null \
            -n {output.sorted}

        echo "== Done =="
        ) > {log} 2>&1
        """
