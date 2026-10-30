# COMEBin: 16 threads, 128 GB RAM, GPU 0; one GPU job at a time until GRES is available.
# Final statuses: success (bins recovered), no_bins (no bins >= 200 kb),
# no_marker_hits (no marker genes detected), no_marker_seeds (no usable marker seeds).
# Other errors, including OOM and timeouts, should fail the job.
rule comebin_sample:
    input:
        contigs="data/assembly/sample/{sample}/{sample}.contigs.fa",
        bam="data/contig_mapping/sample/{sample}/{sample}.sorted.bam",
        bai="data/contig_mapping/sample/{sample}/{sample}.sorted.bam.bai"
    output: status="data/comebin/sample/{sample}/status.tsv"
    params:
        bamdir="data/contig_mapping/sample/{sample}",
        outdir="data/comebin/sample/{sample}",
        seed=42
    threads: 16
    resources: mem_mb=128000, runtime=360, gpu_slot=1, slurm_partition="g1", slurm_extra="--reservation=prj-veo"
    container: config["containers"]["comebin_gpu"]
    log: "logs/binning/comebin/sample/{sample}.log"
    benchmark: "benchmarks/binning/comebin/sample/{sample}.tsv"
    shell:
        """
        if CUDA_VISIBLE_DEVICES=0 run_comebin.sh -a {input.contigs} -p {params.bamdir} -o {params.outdir} -d cuda -s {params.seed} -t {threads} > {log} 2>&1; then
            [[ -f {params.outdir}/comebin_res/comebin_res.tsv && -d {params.outdir}/comebin_res/comebin_res_bins ]] || exit 1
            if [[ -s {params.outdir}/comebin_res/comebin_res.tsv ]]; then status=success; else status=no_bins; fi
        elif grep -q 'HMMsearch completed but found zero marker hits' {log}; then status=no_marker_hits
        elif grep -q 'no usable marker seed set could be generated' {log}; then status=no_marker_seeds
        else exit 1
        fi
        printf 'sample\\tstatus\\n%s\\t%s\\n' '{wildcards.sample}' "$status" > {output.status}
        """

# Prefix sample bin filenames.
rule comebin_sample_rename_bins:
    input: "data/comebin/sample/{sample}/status.tsv"
    output: "data/comebin/sample/{sample}/bins_renamed.done"
    params: bins="data/comebin/sample/{sample}/comebin_res/comebin_res_bins"
    shell:
        r"""
        shopt -s nullglob
        for f in {params.bins:q}/*.fa; do
            name=$(basename "$f")
            [[ "$name" =~ ^[0-9]+\.fa$ ]] || continue
            mv "$f" "$(dirname "$f")/{wildcards.sample}_COMEBin_$name"
        done
        touch {output:q}
        """

# SemiBin2: individual assemblies with multi-sample coverage per biome.
# One A100 GPU, 64 threads, 900 GB RAM, 24 hours; exclusive node allocation.
# Status: success (bins recovered), no_bins (no bins recovered); other errors fail.
rule semibin2_multisample:
    input:
        contigs="data/assembly/multisample/{biome}/{biome}.contigs.fa",
        bams=lambda wc: [f"data/contig_mapping/multisample/{wc.biome}/{s}.sorted.bam" for s in SAMPLES_BY_BIOME[wc.biome]],
        bais=lambda wc: [f"data/contig_mapping/multisample/{wc.biome}/{s}.sorted.bam.bai" for s in SAMPLES_BY_BIOME[wc.biome]]
    output: status="data/semibin2/multisample/{biome}/status.tsv"
    params:
        outdir="data/semibin2/multisample/{biome}",
        samples=lambda wc: SAMPLES_BY_BIOME[wc.biome]
    threads: 64
    resources: mem_mb=900000, runtime=1440, gpu_slot=1, slurm_partition="g1", slurm_extra="--reservation=prj-veo --exclusive"
    container: config["containers"]["semibin2_gpu"]
    log: "logs/binning/semibin2/multisample/{biome}.log"
    benchmark: "benchmarks/binning/semibin2/multisample/{biome}.tsv"
    shell:
        """
        CUDA_VISIBLE_DEVICES=0 SemiBin2 multi_easy_bin -i {input.contigs} -b {input.bams} -o {params.outdir} -s _k --self-supervised --engine gpu --compression none -p {threads} --random-seed 42 > {log} 2>&1
        printf 'sample\\tstatus\\n' > {output.status}
        for s in {params.samples}; do
            [[ -d {params.outdir}/samples/$s/output_bins ]] || exit 1
            if compgen -G "{params.outdir}/samples/$s/output_bins/*.fa" > /dev/null; then status=success; else status=no_bins; fi
            printf '%s\\t%s\\n' "$s" "$status" >> {output.status}
        done
        """

# Restore sample prefixes in multisample bins.
rule semibin2_multisample_restore_headers:
    input: "data/semibin2/multisample/{biome}/status.tsv"
    output: "data/semibin2/multisample/{biome}/headers_restored.done"
    params: bins="data/semibin2/multisample/{biome}/bins"
    shell:
        r"""
        shopt -s nullglob
        for f in {params.bins:q}/*.fa; do
            name=$(basename "$f")
            sample=${{name%%_SemiBin_*}}
            sed -i "/^>/{{ /^>${{sample}}_k/! s/^>/>${{sample}}_k/; }}" "$f"
        done
        touch {output:q}
        """

# VAMB: multi-sample coverage with per-sample binsplitting.
# One A100 GPU, 64 threads, 900 GB RAM, 24 hours.
# Status: success, no_bins; all other errors fail.
rule vamb_multisample:
    input:
        contigs="data/assembly/multisample/{biome}/{biome}.contigs.fa",
        bams=lambda wc: [f"data/contig_mapping/multisample/{wc.biome}/{s}.sorted.bam" for s in SAMPLES_BY_BIOME[wc.biome]],
        bais=lambda wc: [f"data/contig_mapping/multisample/{wc.biome}/{s}.sorted.bam.bai" for s in SAMPLES_BY_BIOME[wc.biome]]
    output:
        status="data/vamb/multisample/{biome}/status.tsv",
        results=directory("data/vamb/multisample/{biome}/run")
    params:
        bamdir="data/contig_mapping/multisample/{biome}",
        samples=lambda wc: SAMPLES_BY_BIOME[wc.biome]
    threads: 64
    resources: mem_mb=900000, runtime=1440, gpu_slot=1, slurm_partition="g1", slurm_extra="--reservation=prj-veo --exclusive"
    container: config["containers"]["vamb_gpu"]
    log: "logs/binning/vamb/multisample/{biome}.log"
    benchmark: "benchmarks/binning/vamb/multisample/{biome}.tsv"
    shell:
        """
        CUDA_VISIBLE_DEVICES=0 vamb bin default --outdir {output.results} --fasta {input.contigs} --bamdir {params.bamdir} -o _k --minfasta 200000 -p {threads} --cuda --seed 42 > {log} 2>&1
        [[ -d {output.results}/bins ]] || exit 1
        printf 'sample\\tstatus\\n' > {output.status}
        for s in {params.samples}; do
            if compgen -G "{output.results}/bins/${{s}}_k*.fna" > /dev/null; then status=success; else status=no_bins; fi
            printf '%s\\t%s\\n' "$s" "$status" >> {output.status}
        done
        """

# Prefix VAMB bin filenames with sample and tool.
rule vamb_multisample_rename_bins:
    input: "data/vamb/multisample/{biome}/run"
    output: "data/vamb/multisample/{biome}/bins_renamed.done"
    shell:
        r"""
        shopt -s nullglob
        for f in {input:q}/bins/*.fna; do
            name=$(basename "$f")
            [[ "$name" =~ ^([0-9]+)_k([0-9]+)\.fna$ ]] || continue
            mv "$f" "$(dirname "$f")/${{BASH_REMATCH[1]}}_VAMB_${{BASH_REMATCH[2]}}.fna"
        done
        touch {output:q}
        """

# SemiBin2: biome coassembly with multi-sample coverage.
rule semibin2_coassembly:
    input:
        contigs="data/assembly/coassembly/{biome}/{biome}.contigs.fa",
        bams=lambda wc: [f"data/contig_mapping/coassembly/{wc.biome}/{s}.sorted.bam" for s in SAMPLES_BY_BIOME[wc.biome]],
        bais=lambda wc: [f"data/contig_mapping/coassembly/{wc.biome}/{s}.sorted.bam.bai" for s in SAMPLES_BY_BIOME[wc.biome]]
    output: bins=directory("data/semibin2/coassembly/{biome}/output_bins")
    params:
        outdir="data/semibin2/coassembly/{biome}",
        seed=42
    threads: 60
    resources: mem_mb=900000, runtime=1440, gpu_slot=1, slurm_partition="g1", slurm_extra="--reservation=prj-veo --exclusive"
    container: config["containers"]["semibin2_gpu"]
    log: "logs/binning/semibin2/coassembly/{biome}.log"
    benchmark: "benchmarks/binning/semibin2/coassembly/{biome}.tsv"
    shell:
        """
        CUDA_VISIBLE_DEVICES=0 SemiBin2 single_easy_bin -i {input.contigs} -b {input.bams} -o {params.outdir} --self-supervised --engine gpu --compression none --minfasta-kbs 200 -p {threads} --random-seed {params.seed} > {log} 2>&1
        """

# Restore biome prefixes in coassembly headers and bin filenames.
rule semibin2_coassembly_restore_headers:
    input: "data/semibin2/coassembly/{biome}/output_bins"
    output: "data/semibin2/coassembly/{biome}/headers_restored.done"
    shell:
        r"""
        shopt -s nullglob
        for f in {input:q}/*.fa; do
            sed -i '/^>k/s/^>/>{wildcards.biome}_/' "$f"
        done
        for f in {input:q}/SemiBin_*.fa; do
            mv "$f" "$(dirname "$f")/{wildcards.biome}_$(basename "$f")"
        done
        touch {output:q}
        """

# COMEBin: biome coassembly with multi-sample coverage.
# Share gpu00 using GPU 1, 48 threads and 550 GB RAM.
rule comebin_coassembly:
    input:
        contigs="data/assembly/coassembly/{biome}/{biome}.contigs.fa",
        bams=lambda wc: [f"data/contig_mapping/coassembly/{wc.biome}/{s}.sorted.bam" for s in SAMPLES_BY_BIOME[wc.biome]],
        bais=lambda wc: [f"data/contig_mapping/coassembly/{wc.biome}/{s}.sorted.bam.bai" for s in SAMPLES_BY_BIOME[wc.biome]]
    output:
        assignments="data/comebin/coassembly/{biome}/comebin_res/comebin_res.tsv",
        bins=directory("data/comebin/coassembly/{biome}/comebin_res/comebin_res_bins")
    params:
        bamdir="data/contig_mapping/coassembly/{biome}",
        outdir="data/comebin/coassembly/{biome}",
        seed=42
    threads: 48
    resources:
        mem_mb=550000,
        runtime=4320,
        qos="long",
        gpu_slot=1,
        slurm_partition="g1",
        slurm_extra="--reservation=prj-veo --nodelist=gpu00"
    container: config["containers"]["comebin_gpu"]
    log: "logs/binning/comebin/coassembly/{biome}.log"
    benchmark: "benchmarks/binning/comebin/coassembly/{biome}.tsv"
    shell:
        """
        CUDA_VISIBLE_DEVICES=1 run_comebin.sh \
            -a {input.contigs} \
            -p {params.bamdir} \
            -o {params.outdir} \
            -d cuda \
            -s {params.seed} \
            -t {threads} \
            > {log} 2>&1
        """

# Prefix coassembly bin filenames.
rule comebin_coassembly_rename_bins:
    input: "data/comebin/coassembly/{biome}/comebin_res/comebin_res_bins"
    output: "data/comebin/coassembly/{biome}/bins_renamed.done"
    shell:
        r"""
        shopt -s nullglob
        for f in {input:q}/*.fa; do
            name=$(basename "$f")
            [[ "$name" =~ ^[0-9]+\.fa$ ]] || continue
            mv "$f" "$(dirname "$f")/{wildcards.biome}_COMEBin_$name"
        done
        touch {output:q}
        """

