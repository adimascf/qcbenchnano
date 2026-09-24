rule align_minimap2:
	input:
		reads=rules.downsample_rasusa.output.reads,
		reference=rules.faidx_mutref.output.fasta,
	log:
		LOGS / "alignment/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.minimap2.log"
	threads: 4
	resources:
		mem="32GiB",
		runtime="1h"
	conda:
		ENVS / "align.yaml"
	params:
		preset="map-ont",
		opts="-aL --cs --MD"
	output:
		bam=temp(RESULTS / "alignment/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.minimap2.sorted.bam"),
		bai=temp(RESULTS / "alignment/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.minimap2.sorted.bam.bai")
	shell:
		"""
		(minimap2 {params.opts} -t {threads} -x {params.preset} {input.reference} {input.reads} | \
				samtools sort -@ {threads} -o {output.bam}) 2> {log}
		samtools index {output.bam} -o {output.bai} 2>> {log}
		"""

rule call_clair3:
	input:
		alignment=rules.align_minimap2.output.bam,
		bai=rules.align_minimap2.output.bai,
		reference=rules.align_minimap2.input.reference,
		faidx=rules.faidx_mutref.output.faidx
	log:
		LOGS / "calling/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.clair3.log"
	threads: 4
	resources:
		mem="128GiB",
		runtime="6h"
	container:
		"docker://hkubal/clair3:v2.0.1"
	shadow:
		"shallow"
	output:
		vcf=RESULTS / "calling/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.clair3.vcf.gz"
	script:
		"../scripts/calling/clair3.sh"
