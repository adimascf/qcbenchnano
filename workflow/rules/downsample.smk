rule downsample_rasusa:
	input:
		reads=RESULTS / "QC/quality/{tool}-{trimmer}/{model}/{sample}.{tool}-{trimmer}.fastq",
		faidx=rules.faidx_self.output.faidx
	log:
		LOGS / "QC/downsampling/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.rasusa.log"
	resources:
		mem="64GiB",
		runtime="30m"
	conda:
		ENVS / "rasusa.yaml"
	params:
		seed="1"
	output:
		reads=RESULTS / "QC/downsampling/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.rasusa.fastq",
	shell:
		"""
		rasusa reads -c {wildcards.depth} -g {input.faidx} -s {params.seed} -o {output.reads} {input.reads} 2> {log}
		"""
