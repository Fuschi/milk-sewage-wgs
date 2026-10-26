#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

usage() {
    cat <<'EOF'
Usage:
  bash quantify_kma_reads.sh \
    --res FILE \
    --mapstat FILE \
    --output-prefix PATH \
    [--min-identity NUMBER] \
    [--min-coverage NUMBER]

Required options:
  --res FILE             KMA results file (.res).
  --mapstat FILE         KMA mapping statistics file (.mapstat).
  --output-prefix PATH   Output path prefix. Parent directory must exist.

Optional options:
  --min-identity NUMBER  Minimum query identity percentage [default: 90].
  --min-coverage NUMBER  Minimum template coverage percentage [default: 60].
  -h, --help             Show this help message.

Outputs:
  PREFIX.all.tsv         All reported references and filter status.
  PREFIX.filtered.tsv    References passing both filters.

Abundance measures:
  aligned_fragments        Aligned fragment count from fragmentCountAln.
  aligned_bases            Aligned base count from per-reference bpTotal.
  fpm                      Fragments per million input fragments.
  fpkm                     Fragments per kilobase per million input fragments.
  tpm_filtered             Length-normalized composition among passing references.
  fragment_proportion_pct  Percentage of input fragments assigned to a reference.
  depth                    Mean reference depth reported by KMA.

TPM is NA for references failing the filters, or when its denominator is zero.
Each row represents one reference sequence; gene variants are not aggregated.
EOF
}

die() {
    echo "Error: $*" >&2
    exit 1
}

require_value() {
    if [[ $# -lt 2 || -z "$2" || "$2" == -* ]]; then
        die "Missing value for $1."
    fi
}

res=""
mapstat=""
prefix=""
min_identity=90
min_coverage=60

while [[ $# -gt 0 ]]; do
    case "$1" in
        --res)
            require_value "$@"
            res=$2
            shift 2
            ;;
        --mapstat)
            require_value "$@"
            mapstat=$2
            shift 2
            ;;
        --output-prefix)
            require_value "$@"
            prefix=$2
            shift 2
            ;;
        --min-identity)
            require_value "$@"
            min_identity=$2
            shift 2
            ;;
        --min-coverage)
            require_value "$@"
            min_coverage=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown argument: $1. Use --help for usage."
            ;;
    esac
done

[[ -n "$res" ]] || die "Missing required option: --res."
[[ -n "$mapstat" ]] || die "Missing required option: --mapstat."
[[ -n "$prefix" ]] || die "Missing required option: --output-prefix."

for file in "$res" "$mapstat"; do
    [[ -f "$file" && -r "$file" && -s "$file" ]] ||
        die "Missing, unreadable, or empty file: $file"
done

for output in "${prefix}.all.tsv" "${prefix}.filtered.tsv"; do
    if [[ "$output" -ef "$res" || "$output" -ef "$mapstat" ]]; then
        die "An output file would overwrite an input file: $output"
    fi
done

awk \
    -v min_identity="$min_identity" \
    -v min_coverage="$min_coverage" \
    -v all_out="${prefix}.all.tsv" \
    -v filtered_out="${prefix}.filtered.tsv" '
BEGIN {
    FS = OFS = "\t"

    if (min_identity !~ /^[0-9]+([.][0-9]+)?$/ ||
        min_coverage !~ /^[0-9]+([.][0-9]+)?$/ ||
        min_identity > 100 || min_coverage > 100) {
        fail("Thresholds must be numbers between 0 and 100.")
    }
}

function fail(message) {
    failed = 1
    print "Error: " message > "/dev/stderr"
    exit 1
}

function clean(value) {
    sub(/^[[:space:]]+/, "", value)
    sub(/[[:space:]]+$/, "", value)
    return value
}

function is_number(value) {
    return clean(value) ~ /^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$/
}

# Read metadata and reference statistics from the mapstat file.
FNR == NR {
    if ($1 == "## fragmentCount") {
        if (!is_number($2)) {
            fail("Invalid total fragment count.")
        }

        total_fragments = $2 + 0
        next
    }

    if ($1 ~ /^# refSequence/) {
        for (i = 1; i <= NF; i++) {
            map_col[clean($i)] = i
        }

        if (!("fragmentCountAln" in map_col) ||
            !("bpTotal" in map_col)) {
            fail("Missing fragmentCountAln or bpTotal column.")
        }

        map_header_found = 1
        next
    }

    if ($0 ~ /^#/ || $0 ~ /^[[:space:]]*$/) {
        next
    }

    if (!map_header_found) {
        fail("Unrecognized mapstat header.")
    }

    reference = clean($1)

    if (reference in fragments) {
        fail("Duplicate mapstat reference: " reference)
    }

    if (!is_number($(map_col["fragmentCountAln"])) ||
        !is_number($(map_col["bpTotal"]))) {
        fail("Invalid mapping statistics for: " reference)
    }

    fragments[reference] = $(map_col["fragmentCountAln"]) + 0
    aligned_bases[reference] = $(map_col["bpTotal"]) + 0
    next
}

# Read and validate the results header.
FNR == 1 {
    if (!map_header_found || total_fragments <= 0) {
        fail("Missing mapstat metadata or invalid total fragment count.")
    }

    for (i = 1; i <= NF; i++) {
        name = clean($i)
        sub(/^#/, "", name)
        res_col[name] = i
    }

    required = "Template Template_length Template_Identity Template_Coverage Query_Identity Depth"
    n_required = split(required, names, " ")

    for (i = 1; i <= n_required; i++) {
        if (!(names[i] in res_col)) {
            fail("Missing results column: " names[i])
        }
    }

    res_header_found = 1
    next
}

# Store reference results and calculate the filtered TPM denominator.
$0 !~ /^#/ && $0 !~ /^[[:space:]]*$/ {
    reference = clean($(res_col["Template"]))

    if (!(reference in fragments)) {
        fail("Reference missing from mapstat: " reference)
    }

    if (reference in seen) {
        fail("Duplicate results reference: " reference)
    }

    seen[reference] = 1

    for (i = 2; i <= n_required; i++) {
        if (!is_number($(res_col[names[i]]))) {
            fail("Invalid " names[i] " value for: " reference)
        }
    }

    length_bp = $(res_col["Template_length"]) + 0

    if (length_bp <= 0) {
        fail("Invalid reference length: " reference)
    }

    n++
    references[n] = reference
    lengths[n] = length_bp
    template_identity[n] = $(res_col["Template_Identity"]) + 0
    query_identity[n] = $(res_col["Query_Identity"]) + 0
    coverage[n] = $(res_col["Template_Coverage"]) + 0
    depth[n] = $(res_col["Depth"]) + 0

    passed[n] = (query_identity[n] >= min_identity && coverage[n] >= min_coverage)

    if (passed[n]) {
        tpm_denominator += fragments[reference] / length_bp
        n_pass++
    }
}

# Write outputs after all input records have been validated.
END {
    if (failed) {
        exit 1
    }

    if (!res_header_found) {
        print "Error: Missing results header." > "/dev/stderr"
        exit 1
    }

    header = "reference" OFS "length_bp" OFS \
             "template_identity_pct" OFS "query_identity_pct" OFS \
             "template_coverage_pct" OFS "depth" OFS \
             "aligned_fragments" OFS "aligned_bases" OFS \
             "total_fragments" OFS "fpm" OFS "fpkm" OFS \
             "tpm_filtered" OFS "fragment_proportion_pct" OFS \
             "pass_filters"

    print header > all_out
    print header > filtered_out

    for (i = 1; i <= n; i++) {
        reference = references[i]
        count = fragments[reference]

        fpm = count * 1000000 / total_fragments
        fpkm = count * 1000000000 / (lengths[i] * total_fragments)
        proportion = count * 100 / total_fragments

        # TPM uses only references passing both filters.
        tpm = "NA"

        if (passed[i] && tpm_denominator > 0) {
            tpm = sprintf("%.6f", (count / lengths[i]) * 1000000 / tpm_denominator)
        }

        row = reference OFS lengths[i] OFS \
              template_identity[i] OFS query_identity[i] OFS \
              coverage[i] OFS depth[i] OFS \
              sprintf("%.0f", count) OFS \
              sprintf("%.0f", aligned_bases[reference]) OFS \
              sprintf("%.0f", total_fragments) OFS \
              sprintf("%.6f", fpm) OFS \
              sprintf("%.6f", fpkm) OFS \
              tpm OFS sprintf("%.8f", proportion) OFS passed[i]

        print row > all_out

        if (passed[i]) {
            print row > filtered_out
        }
    }

    if (close(all_out) != 0) {
        fail("Failed to write output: " all_out)
    }

    if (close(filtered_out) != 0) {
        fail("Failed to write output: " filtered_out)
    }

    print "References processed: " n + 0 > "/dev/stderr"
    print "References passing filters: " n_pass + 0 > "/dev/stderr"
}
' "$mapstat" "$res"
