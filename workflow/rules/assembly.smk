REPEAT = int(config.get("repeat", 1))

rule assembly_flye:
	input:
		reads=rules.downsample_rasusa.output.reads
	log:
		LOGS / "assembly/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.flye.log"
	threads: 16
	resources:
		mem="128GiB",
		runtime=f"{10 * REPEAT}h"
	conda:
		ENVS / "flye.yaml"
	params:
		ont="--nano-hq",
	output:
		assembly=RESULTS / "assembly/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.{depth}x.assembly.fasta",
		graph=RESULTS / "assembly/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.{depth}x.assembly_graph.gfa",
		info=RESULTS / "assembly/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.{depth}x.assembly_info.txt"
	benchmark:
		repeat(BENCHMARK / "assembly/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.tsv", REPEAT)
	shell:
		"""
		tmp_results=$(mktemp -d)
		flye --debug {params.ont} {input.reads} --out-dir $tmp_results --threads {threads} 2> {log}
		mv "${{tmp_results}}"/assembly.fasta {output.assembly} 2>> {log}
		mv "${{tmp_results}}"/assembly_graph.gfa {output.graph} 2>> {log}
		mv "${{tmp_results}}"/assembly_info.txt {output.info} 2>> {log}

		rm -rf $tmp_results
		"""
	
rule reorient_assembly_sample:
	input:
		assembly=rules.assembly_flye.output.assembly
	log:
		LOGS / "assembly/{tool}-{trimmer}/{depth}x/{model}/{sample}.{tool}-{trimmer}.dnaapler.log"
	threads: 8
	resources:
		mem="64GiB",
		runtime=f"{15* REPEAT}h"
	conda:
		ENVS / "dnaapler.yaml"
	params:
		seed="--seed_value 8",
		prefix="{sample}.{tool}-{trimmer}.{depth}x.assembly",
	output:
		assembly=RESULTS / "assembly/{tool}-{trimmer}/{depth}x/{model}/dnaapler/{sample}.{tool}-{trimmer}.{depth}x.assembly_reoriented.fasta"
	shell:
		"""
		# create a unique temporary directory
		tmp_results=$(mktemp -d)
		
		# Run dnaapler inside the unique temp directory
		dnaapler all {params.seed} -f -i {input.assembly} -o $tmp_results -p {params.prefix} -t {threads} 2> {log}
		
		# Move the final fasta to your actual Snakemake output target
		mv "${{tmp_results}}/{params.prefix}_reoriented.fasta" {output.assembly} 2>> {log}
		
		# Clean up the temp directory
		rm -rf $tmp_results
		"""

use rule reorient_assembly_sample as reorient_assembly_reference with:
	input:
		assembly=get_reference_genome
	log:
		LOGS / "reference/{sample}_dnaapler.log"
	params:
		prefix="{sample}",
	output:
		assembly=RESULTS / "reference/dnaapler/{sample}_reoriented.fasta"

