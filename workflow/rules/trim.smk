REPEAT = int(config.get("repeat", 1))

# trimmer = "fastplong"
# rule trim_fastplong:
# 	input:
# 		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
# 	log:
# 		LOGS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.log"
# 	threads: 8
# 	resources:
# 		mem="64GiB",
# 		runtime=f"{4 * REPEAT}m",
# 	conda:
# 		ENVS / "fastplong.yaml"
# 	params:
# 		nofilter="--disable_quality_filtering --disable_length_filtering"
# 	output:
# 		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.fastq"),
# 		json=temp(RESULTS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.json"),
# 		html=temp(RESULTS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.html")
# 	benchmark:
# 		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.tsv", REPEAT)
# 	shell:
# 		"""
# 		fastplong -i {input.reads} -o {output.reads} {params.nofilter} --thread {threads} \
# 				--json {output.json} --html {output.html} --verbose 2> {log}
# 		"""

# trimmer = "porechop_abi_split"
# # by default, it will split reads when an adapter is found in the middle
# rule trim_porcechop_abi_split:
# 	input:
# 		reads=lambda wildcards: get_original_fastqs(wildcards, "trim") 
# 	log: 
# 		LOGS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.log"
# 	threads: 32 
# 	resources:
# 		mem="256GiB",
# 		runtime=f"{24 * REPEAT}h"
# 	conda:
# 		ENVS / "porechop_abi.yaml"
# 	output:
# 		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.fastq")
# 	benchmark: 
# 		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.tsv", REPEAT)
# 	shell: 
# 		"""
# 		porechop_abi -abi -t {threads} -i {input.reads} -o {output.reads} 2> {log}
# 		"""
#
trimmer = "porechop_abi" # "porechop_abi_discard"
rule trim_porcechop_abi_discard:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
	log:
		LOGS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.log"
	threads: 32
	resources:
		mem="256GiB",
		runtime=f"{24 * REPEAT}h"
	conda:
		ENVS / "porechop_abi.yaml"
	params:
		nochimera="--discard_middle"
	output:
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.fastq")
	benchmark: 
		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.tsv", REPEAT)
	shell: 
		"""
		porechop_abi -abi -t {threads} {params.nochimera} -i {input.reads} -o {output.reads} 2> {log}
		"""

# trimmer = "porechop_abi_nocheck"
# rule trim_porcechop_abi_nocheck:
# 	input:
# 		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
# 	log:
# 		LOGS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.log"
# 	threads: 32
# 	resources:
# 		mem="256GiB",
# 		runtime=f"{24 * REPEAT}h"
# 	conda:
# 		ENVS / "porechop_abi.yaml"
# 	params:
# 		nocheck="--no_split" # will skip chimera searching
# 	output:
# 		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.fastq")
# 	benchmark: 
# 		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.tsv", REPEAT)
# 	shell:  
# 		""" 
# 		porechop_abi -abi -t {threads} {params.nocheck} -i {input.reads} -o {output.reads} 2> {log}
# 		"""  

trimmer = "dorado"
rule trim_dorado:  
	input:  
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")  
	log:  
		LOGS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.log"  
	resources:  
		mem="32GiB",  
		runtime=f"{15 * REPEAT}m"  
	container:  
		"docker://nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"
	params:  
		kit=lambda wildcards: get_sequencing_kits(wildcards),  
		output_fq="--emit-fastq"  
	output:  
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.fastq")
	benchmark:  
		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.tsv", REPEAT)
	shell: 
		"""
		(dorado trim {params.output_fq} --sequencing-kit {params.kit} {input.reads} > {output.reads}) 2>> {log}
		"""

trimmer = "barbell" #"barbell_default"
rule trim_barbell:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "notrim")
	log:
		LOGS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.log"
	resources:
		mem="32GiB",
		runtime=f"{20 * REPEAT}h"
	conda:
		ENVS / "barbell_seqkit.yaml"
	params:
		kit=lambda wildcards: get_sequencing_kits(wildcards),
		maximize="no"
	output:
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.tsv", REPEAT)
	script:
		"../scripts/trimming/barbell_kit.sh"

# trimmer = "barbell_max"
# rule barbell_trim_max:
# 	input:
# 		reads=lambda wildcards: get_original_fastqs(wildcards, "notrim")
# 	log:
# 		LOGS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.log"
# 	resources:
# 		mem="32GiB",
# 		runtime=f"{20 * REPEAT}h"
# 	conda:
# 		ENVS / "barbell_seqkit.yaml"
# 	params:
# 		kit=lambda wildcards: get_sequencing_kits(wildcards),
# 		maximize="yes"
# 	output:
# 		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.fastq")
# 	benchmark:
# 		repeat(BENCHMARK / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.tsv", REPEAT)
# 	script:
# 		"../scripts/trimming/barbell_kit.sh"

trimmer = "untrimmed"
rule trim_notrim:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "notrim")
	log:
		LOGS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.log"
	resources:
		mem="8GiB",
		runtime="10m"
	output:
		reads=temp(RESULTS / f"QC/trimming/{trimmer}/{{model}}/{{sample}}.{trimmer}.fastq")
	shell:
		"""
		# here we just copy the input file, no trimming process being perfomed
		cp {input.reads} {output.reads}
		"""
