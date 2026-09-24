REPEAT = int(config.get("repeat", 1))

# 1. Fastplong (adapter trimming with quality/length filtering disabled)
rule trim_fastplong:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
	log:
		LOGS / "QC/trimming/fastplong/{model}/{sample}.log"
	threads: 8
	resources:
		mem="64GiB",
		runtime=f"{4 * REPEAT}m"
	conda:
		ENVS / "fastplong.yaml"
	params:
		nofilter="--disable_quality_filtering --disable_length_filtering"
	output:
		reads=temp(RESULTS / "QC/trimming/fastplong/{model}/{sample}.fastplong.fastq"),
		json=temp(RESULTS / "QC/trimming/fastplong/{model}/{sample}.fastplong.json"),
		html=temp(RESULTS / "QC/trimming/fastplong/{model}/{sample}.fastplong.html")
	benchmark:
		repeat(BENCHMARK / "QC/trimming/fastplong/{model}/{sample}.fastplong.tsv", REPEAT)
	shell:
		"""
		fastplong -i {input.reads} -o {output.reads} {params.nofilter} --thread {threads} \
				--json {output.json} --html {output.html} --verbose 2> {log}
		"""

# 2. Porechop_ABI split (default splitting mode)
rule trim_porechop_abi_split:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
	log:
		LOGS / "QC/trimming/porechop_abi_split/{model}/{sample}.log"
	threads: 32
	resources:
		mem="256GiB",
		runtime=f"{24 * REPEAT}h"
	conda:
		ENVS / "porechop_abi.yaml"
	output:
		reads=temp(RESULTS / "QC/trimming/porechop_abi_split/{model}/{sample}.porechop_abi_split.fastq")
	benchmark:
		repeat(BENCHMARK / "QC/trimming/porechop_abi_split/{model}/{sample}.porechop_abi_split.tsv", REPEAT)
	shell:
		"""
		porechop_abi -abi -t {threads} -i {input.reads} -o {output.reads} 2> {log}
		"""

# 3. Porechop_ABI discard (--discard_middle)
rule trim_porechop_abi_discard:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
	log:
		LOGS / "QC/trimming/porechop_abi_discard/{model}/{sample}.log"
	threads: 32
	resources:
		mem="256GiB",
		runtime=f"{24 * REPEAT}h"
	conda:
		ENVS / "porechop_abi.yaml"
	params:
		nochimera="--discard_middle"
	output:
		reads=temp(RESULTS / "QC/trimming/porechop_abi_discard/{model}/{sample}.porechop_abi_discard.fastq")
	benchmark:
		repeat(BENCHMARK / "QC/trimming/porechop_abi_discard/{model}/{sample}.porechop_abi_discard.tsv", REPEAT)
	shell:
		"""
		porechop_abi -abi -t {threads} {params.nochimera} -i {input.reads} -o {output.reads} 2> {log}
		"""

# 4. Porechop_ABI nocheck (--no_split)
rule trim_porechop_abi_nocheck:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
	log:
		LOGS / "QC/trimming/porechop_abi_nocheck/{model}/{sample}.log"
	threads: 32
	resources:
		mem="256GiB",
		runtime=f"{24 * REPEAT}h"
	conda:
		ENVS / "porechop_abi.yaml"
	params:
		nocheck="--no_split"
	output:
		reads=temp(RESULTS / "QC/trimming/porechop_abi_nocheck/{model}/{sample}.porechop_abi_nocheck.fastq")
	benchmark:
		repeat(BENCHMARK / "QC/trimming/porechop_abi_nocheck/{model}/{sample}.porechop_abi_nocheck.tsv", REPEAT)
	shell:
		"""
		porechop_abi -abi -t {threads} {params.nocheck} -i {input.reads} -o {output.reads} 2> {log}
		"""

# 5. Dorado (trimmed during basecalling)
rule trim_dorado:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "trim")
	log:
		LOGS / "QC/trimming/dorado/{model}/{sample}.log"
	resources:
		mem="32GiB",
		runtime=f"{15 * REPEAT}m"
	container:
		"docker://nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"
	params:
		kit=lambda wildcards: get_sequencing_kits(wildcards),
		output_fq="--emit-fastq"
	output:
		reads=temp(RESULTS / "QC/trimming/dorado/{model}/{sample}.dorado.fastq")
	benchmark:
		repeat(BENCHMARK / "QC/trimming/dorado/{model}/{sample}.dorado.tsv", REPEAT)
	shell:
		"""
		(dorado trim {params.output_fq} --sequencing-kit {params.kit} {input.reads} > {output.reads}) 2>> {log}
		"""

# 6. Barbell default
rule trim_barbell_default:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "notrim")
	log:
		LOGS / "QC/trimming/barbell_default/{model}/{sample}.log"
	resources:
		mem="32GiB",
		runtime=f"{20 * REPEAT}h"
	conda:
		ENVS / "barbell_seqkit.yaml"
	params:
		kit=lambda wildcards: get_sequencing_kits(wildcards),
		maximize="no"
	output:
		reads=temp(RESULTS / "QC/trimming/barbell_default/{model}/{sample}.barbell_default.fastq")
	benchmark:
		repeat(BENCHMARK / "QC/trimming/barbell_default/{model}/{sample}.barbell_default.tsv", REPEAT)
	script:
		"../scripts/trimming/barbell_kit.sh"

# 7. Barbell max (--maximize)
rule trim_barbell_max:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "notrim")
	log:
		LOGS / "QC/trimming/barbell_max/{model}/{sample}.log"
	resources:
		mem="32GiB",
		runtime=f"{20 * REPEAT}h"
	conda:
		ENVS / "barbell_seqkit.yaml"
	params:
		kit=lambda wildcards: get_sequencing_kits(wildcards),
		maximize="yes"
	output:
		reads=temp(RESULTS / "QC/trimming/barbell_max/{model}/{sample}.barbell_max.fastq")
	benchmark:
		repeat(BENCHMARK / "QC/trimming/barbell_max/{model}/{sample}.barbell_max.tsv", REPEAT)
	script:
		"../scripts/trimming/barbell_kit.sh"

# 8. Untrimmed
rule trim_notrim:
	input:
		reads=lambda wildcards: get_original_fastqs(wildcards, "notrim")
	log:
		LOGS / "QC/trimming/untrimmed/{model}/{sample}.log"
	resources:
		mem="8GiB",
		runtime="10m"
	output:
		reads=temp(RESULTS / "QC/trimming/untrimmed/{model}/{sample}.untrimmed.fastq")
	shell:
		"""
		cp {input.reads} {output.reads}
		"""
