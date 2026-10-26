# COMEBin: sample assignments.
rule dastool_comebin_sample_input:
    input: "data/binning/comebin/sample/{sample}/bins_renamed.done"
    output: "data/refinement/dastool/sample/{sample}/inputs/comebin.tsv"
    params: bins="data/binning/comebin/sample/{sample}/comebin_res/comebin_res_bins"
    container: config["containers"]["dastool"]
    shell:
        """
        if compgen -G "{params.bins}/*.fa" > /dev/null; then
            Fasta_to_Contig2Bin.sh -i {params.bins} -e fa > {output}
        else
            : > {output}
        fi
        """

# SemiBin2 and VAMB: sample assignments from pooled bins.
rule dastool_multisample_inputs:
    input:
        semibin2=lambda wc: f"data/binning/semibin2/multisample/{BIOME_BY_SAMPLE[wc.sample]}/headers_restored.done",
        vamb=lambda wc: f"data/binning/vamb/multisample/{BIOME_BY_SAMPLE[wc.sample]}/bins_renamed.done"
    output:
        semibin2="data/refinement/dastool/sample/{sample}/inputs/semibin2.tsv",
        vamb="data/refinement/dastool/sample/{sample}/inputs/vamb.tsv"
    params:
        semibin2=lambda wc: f"data/binning/semibin2/multisample/{BIOME_BY_SAMPLE[wc.sample]}/bins",
        vamb=lambda wc: f"data/binning/vamb/multisample/{BIOME_BY_SAMPLE[wc.sample]}/run/bins"
    container: config["containers"]["dastool"]
    shell:
        """
        Fasta_to_Contig2Bin.sh -i {params.semibin2} -e fa | awk -v p="{wildcards.sample}_" 'index($1, p) == 1' > {output.semibin2}
        Fasta_to_Contig2Bin.sh -i {params.vamb} -e fna | awk -v p="{wildcards.sample}_" 'index($1, p) == 1' > {output.vamb}
        """

# COMEBin and SemiBin2: coassembly assignments.
rule dastool_coassembly_inputs:
    input:
        comebin="data/binning/comebin/coassembly/{biome}/bins_renamed.done",
        semibin2="data/binning/semibin2/coassembly/{biome}/headers_restored.done"
    output:
        comebin="data/refinement/dastool/coassembly/{biome}/inputs/comebin.tsv",
        semibin2="data/refinement/dastool/coassembly/{biome}/inputs/semibin2.tsv"
    params:
        comebin="data/binning/comebin/coassembly/{biome}/comebin_res/comebin_res_bins",
        semibin2="data/binning/semibin2/coassembly/{biome}/output_bins"
    container: config["containers"]["dastool"]
    shell:
        """
        Fasta_to_Contig2Bin.sh -i {params.comebin} -e fa > {output.comebin}
        Fasta_to_Contig2Bin.sh -i {params.semibin2} -e fa > {output.semibin2}
        """

# Refine sample bins from COMEBin, SemiBin2 and VAMB.
# Refine sample bins and record runs without results.
rule dastool_sample:
    input:
        contigs="data/assembly/sample/{sample}/{sample}.contigs.fa",
        comebin="data/refinement/dastool/sample/{sample}/inputs/comebin.tsv",
        semibin2="data/refinement/dastool/sample/{sample}/inputs/semibin2.tsv",
        vamb="data/refinement/dastool/sample/{sample}/inputs/vamb.tsv"
    output:
        summary="data/refinement/dastool/sample/{sample}/{sample}_DASTool_summary.tsv",
        assignments="data/refinement/dastool/sample/{sample}/{sample}_DASTool_contig2bin.tsv",
        bins=directory("data/refinement/dastool/sample/{sample}/{sample}_DASTool_bins"),
        status="data/refinement/dastool/sample/{sample}/status.tsv"
    params: prefix="data/refinement/dastool/sample/{sample}/{sample}"
    threads: 16
    resources: mem_mb=32000, runtime=240
    container: config["containers"]["dastool"]
    log: "logs/refinement/dastool/sample/{sample}.log"
    benchmark: "benchmarks/refinement/dastool/sample/{sample}.tsv"
    shell:
        """
        if [[ ! -s {input.comebin} && ! -s {input.semibin2} && ! -s {input.vamb} ]]; then
            status=no_input_bins
            echo "No input bins available." > {log}
        elif DAS_Tool -i {input.comebin},{input.semibin2},{input.vamb} -l comebin,semibin2,vamb -c {input.contigs} -o {params.prefix} -t {threads} --search_engine diamond --score_threshold 0.5 --write_bins --write_bin_evals > {log} 2>&1; then
            status=success
        elif grep -Fq "No bins with bin-score" {log}; then
            status=no_bins_above_threshold
        else
            exit 1
        fi
        if [[ "$status" != success ]]; then
            : > {output.summary}
            : > {output.assignments}
            mkdir -p {output.bins}
        fi
        printf 'sample\\tstatus\\n%s\\t%s\\n' '{wildcards.sample}' "$status" > {output.status}
        """

# Refine coassembly bins from COMEBin and SemiBin2.
rule dastool_coassembly:
    input:
        contigs="data/assembly/coassembly/{biome}/{biome}.contigs.fa",
        comebin="data/refinement/dastool/coassembly/{biome}/inputs/comebin.tsv",
        semibin2="data/refinement/dastool/coassembly/{biome}/inputs/semibin2.tsv"
    output:
        summary="data/refinement/dastool/coassembly/{biome}/{biome}_DASTool_summary.tsv",
        assignments="data/refinement/dastool/coassembly/{biome}/{biome}_DASTool_contig2bin.tsv",
        bins=directory("data/refinement/dastool/coassembly/{biome}/{biome}_DASTool_bins")
    params: prefix="data/refinement/dastool/coassembly/{biome}/{biome}"
    threads: 32
    resources: mem_mb=128000, runtime=1440
    container: config["containers"]["dastool"]
    log: "logs/refinement/dastool/coassembly/{biome}.log"
    benchmark: "benchmarks/refinement/dastool/coassembly/{biome}.tsv"
    shell:
        """
        DAS_Tool -i {input.comebin},{input.semibin2} -l comebin,semibin2 -c {input.contigs} -o {params.prefix} -t {threads} --search_engine diamond --score_threshold 0.5 --write_bins --write_bin_evals > {log} 2>&1
        """
