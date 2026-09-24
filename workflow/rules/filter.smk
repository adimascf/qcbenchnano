REPEAT = config.get("repeat", 1)

rule filter_fastplong:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{sample}.fastq"
	log:
		LOGS / "QC/filtering/fastplong/{trimmer}/{sample}.{trimmer}.fastplong.log"
	threads: 8
	resources:
		mem_mb=64000,
		runtime=f"{5 * REPEAT}m"
	conda: 
		ENVS / "fastplong.yaml"
	params:
		notrim="--disable_quality_filtering" 
	output:
		reads=RESULTS / "QC/filtering/fastplong/{trimmer}/{sample}.{trimmer}.fastplong.fastq",
		json=RESULTS / "QC/filtering/fastplong/{trimmer}/{sample}.{trimmer}.fastplong.json",
		html=RESULTS / "QC/filtering/fastplong/{trimmer}/{sample}.{trimmer}.fastplong.html"
	benchmark:
		repeat(BENCHMARK / "QC/filtering/fastplong/{trimmer}/{sample}.{trimmer}.fastplong.tsv", REPEAT)
	shell:
		"""
		fastplong -i {input.reads} -o {output.reads} {params.notrim} --thread {threads} \
				--json {output.json} --html {output.html} 2> {log}
		"""

rule filter_filtlong:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{sample}.fastq"
	log:
		LOGS / "QC/filtering/fastplong/{trimmer}/{sample}.{trimmer}.filtlong.log"
	resources:
		mem_mb=128000,
		runtime=f"{20 * REPEAT}m"
	conda:
		ENVS / "filtlong.yaml"
	params:
		default="--min_length 1kb --keep_percent 90 --target_bases 500m"
	output:
		reads=RESULTS / "QC/filtering/filtlong/{trimmer}/{sample}.{trimmer}.filtlong.fastq",
	benchmark:
		repeat(BENCHMARK / "QC/filtering/filtlong/{trimmer}/{sample}.{trimmer}.filtlong.tsv", REPEAT)
	shell:
		"""
		filtlong {params.default} {input.reads} | gzip > {output.reads}
		"""

