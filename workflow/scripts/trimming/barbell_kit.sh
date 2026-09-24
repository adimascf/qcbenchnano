#!/usr/bin/env bash
set -euo pipefail

# diret the stderr to the log file
exec 2> "${snakemake_log[0]}"

reads="${snakemake_input[reads]}"
kit="${snakemake_params[kit]}"
outdir=$(mktemp -d)
output="${snakemake_output[reads]}"
max="${snakemake_params[maximize]}"

if [[ "$max" == "yes" ]]; then
	 
	barbell kit --maximize --kit "$kit" --input "$reads" --output "$outdir"

	# concatenate all of the "demultiplexed" output files, and the we remove the duplicate by read id
	cat "$outdir"/*.fastq | seqkit rmdup -o "$output"

elif [[ "$max" == "no" ]]; then
 
	barbell kit --kit "$kit" --input "$reads" --output "$outdir"

	# concatenate all of the "demultiplexed" output files, and the we remove the duplicate by read id
	cat "$outdir"/*.fastq | seqkit rmdup -o "$output"
fi

