#!/usr/bin/env bash
set -euo pipefail

files=$1
outdir=$2

mkdir -p "$outdir"

for file in `ls "$files"/*bam`; do
	base=$(basename "$file")
	name=$(echo "$base" | awk -v FS="." '{print $1}')
	echo "converting $name..."
	samtools fastq  -@ 4 "$file" > "$outdir/$name.fastq"
	rm "$file"
done
