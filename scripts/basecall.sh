#!/usr/bin/env bash
set -euo pipefail


dir=$1
model=$2
outdir=$3
trim_flag=${4:-""} # Default is empty (trimming enabled)

case "$model" in
    "hac")
        model_version="dna_r10.4.1_e8.2_400bps_hac@v5.2.0"
        ;;
    "sup")
        model_version="dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
        ;;
    *)
		echo "Error: Model must be 'hac' or 'sup'. You provided '$model'." >&2
        exit 1
        ;;
esac

# Add log so that the paramters are captured
echo "================================================================="
echo "DORADO BASECALLING START"
echo "Input Directory : $dir"
echo "Explicit Model  : $model_version"
echo "Output Directory: $outdir"
echo "Trim Flag       : ${trim_flag:-"Enabled (Default)"}"
echo "================================================================="


mkdir -p "$outdir"

for run in "$dir"/*/; do
    name=$(basename "${run%/}")
	kit_flag=""

	# Only select a kit if we are NOT skipping trimming
    if [[ "$trim_flag" != "--no-trim" ]]; then
        case "$name" in
			"ATCC_10708__202309"|"ATCC_17802__202309"|"ATCC_19119__202309"|"ATCC_25922__202309"|"ATCC_35897__202309"|"ATCC_BAA-679__202309")
                kit_flag="--kit-name SQK-RBK114-96"
                ;;
            "ATCC_33560__202309"|"ATCC_35221__202309"|"AJ292__202310"|"RDH275__202311"|"MMC234__202311"|"BPH2947__202310"|"AMtb_1__202402")
                kit_flag="--kit-name SQK-NBD114-96"
                ;;
            *)
                echo "Warning: No kit matched for $name. Running without kit flag."
                ;;
        esac
    fi

    dorado basecaller $trim_flag $kit_flag --recursive "$model_version" "$run" > "$outdir/$name.bam"
done
# ssubmit -t 1d -m 16g basecall "./scripts/basecall.sh <inputdir> <model>" -- --qos gpu -p gpu_cuda --gres=gpu:h100:1
