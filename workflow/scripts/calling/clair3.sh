#!/usr/bin/env bash
set -euo pipefail

# send all stderr and stdout from this script to the log file
exec &>"${snakemake_log[0]}"

aln="${snakemake_input[alignment]}"
ref="${snakemake_input[reference]}"
outvcf="${snakemake_output[vcf]}"

sample="${snakemake_wildcards[sample]}"
basecall_model="${snakemake_wildcards[model]}"

if [[ "$basecall_model" == "sup" ]]; then
    model_name="r1041_e82_400bps_sup_v520"
elif [[ "$basecall_model" == "hac" ]]; then
    model_name="r1041_e82_400bps_hac_v520"
else
    echo "Error: Unknown basecalling model '${basecall_model}' provided by Snakemake." >&2
    exit 1
fi

model_path="/scratch/user/s4897040/biox7021/data/clair3_models/${model_name}"
tmpoutdir=$(mktemp -d)

trap 'rm -rf "$tmpoutdir"' EXIT

/opt/bin/run_clair3.sh \
    --bam_fn="$aln" \
    --ref_fn="$ref" \
    --threads="${snakemake[threads]}" \
    --platform="ont" \
    --model_path="$model_path" \
    --output="$tmpoutdir" \
    --sample_name="$sample" \
    --include_all_ctgs \
    --haploid_precise \
    --no_phasing_for_fa \
    --enable_long_indel

mv "${tmpoutdir}/merge_output.vcf.gz" "$outvcf"
