REPEAT = int(config.get("repeat", 1))

rule faidx_self:
	input:
		reference=get_reference_genome,
	log:
		LOGS / "faidx_reference/{sample}.log",
	resources:
		mem="1GiB",
		runtime="5m"
	conda:
		ENVS / "rasusa.yaml"
	output:
		fasta=RESULTS / "reference/{sample}.fa",
		faidx=RESULTS / "reference/{sample}.fa.fai",
	shell:
		"""
		cp {input.reference} {output.fasta} 2> {log}
		samtools faidx {output.fasta} 2>> {log}
		"""

use rule faidx_self as faidx_mutref with: 
	input:
		reference=get_mutreference_genome
	log:
		LOGS / "faidx_reference/{sample}.mutref.log"
	output:
		fasta=RESULTS / "mutreference/{sample}.mutref.fa",
		faidx=RESULTS / "mutreference/{sample}.mutref.fa.fai"

# 1. Fastplong 100 bp (adapter trimming disabled)
tool = "fastplong_100"
rule quality_fastplong_100:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq"
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="16GiB",
		runtime=f"{5 * REPEAT}m"
	conda:
		ENVS / "fastplong.yaml"
	params:
		no_adapt_trimming="--disable_adapter_trimming",
		quality_threshold=get_quality_threshold, # hac is 10, sup is 15
		length_threshold=100
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq"),
		json=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.json"),
		html=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.html")
	benchmark:   
		repeat(BENCHMARK / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.tsv", REPEAT)
	shell:    
		"""   
		fastplong -i {input.reads} -o {output.reads} {params.no_adapt_trimming} \
				--length_required {params.length_threshold} \
				--mean_qual {params.quality_threshold} \
				--json {output.json} --html {output.html} --verbose 2> {log} 
		"""

# Fastplong 1000 bp (adapter trimming disabled)
tool = "fastplong_1000"
rule quality_fastplong_1000:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq"
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="16GiB",
		runtime=f"{5 * REPEAT}m"
	conda:
		ENVS / "fastplong.yaml"
	params:
		no_adapt_trimming="--disable_adapter_trimming",
		quality_threshold=get_quality_threshold, # hac is 10, sup is 15
		length_threshold=1000
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq"),
		json=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.json"),
		html=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.html")
	benchmark:   
		repeat(BENCHMARK / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.tsv", REPEAT)
	shell:    
		"""   
		fastplong -i {input.reads} -o {output.reads} {params.no_adapt_trimming} \
				--length_required {params.length_threshold} \
				--mean_qual {params.quality_threshold} \
				--json {output.json} --html {output.html} --verbose 2> {log} 
		"""

# 2. SeqKit 100 bp
tool = "seqkit_100"
rule quality_seqkit_100:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq"
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="16GiB",
		runtime=f"{5 * REPEAT}m"
	conda:
		ENVS / "seqkit.yaml"
	params:
		quality_threshold=get_quality_threshold, # hac is 10, sup is 15
		length_threshold=100
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq")
	benchmark:   
		repeat(BENCHMARK / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.tsv", REPEAT)
	shell:    
		"""
		seqkit seq --min-len {params.length_threshold} --min-qual {params.quality_threshold} -o {output.reads} {input.reads} 2> {log}
		"""

# SeqKit 1000 bp
tool = "seqkit_1000"
rule quality_seqkit_1000:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq"
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="16GiB",
		runtime=f"{5 * REPEAT}m"
	conda:
		ENVS / "seqkit.yaml"
	params:
		quality_threshold=get_quality_threshold, # hac is 10, sup is 15
		length_threshold=1000
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq")
	benchmark:   
		repeat(BENCHMARK / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.tsv", REPEAT)
	shell:    
		"""
		seqkit seq --min-len {params.length_threshold} --min-qual {params.quality_threshold} -o {output.reads} {input.reads} 2> {log}
		"""

# 3. Filtlong Mean Quality Only (weight 0 for length and window_q)
tool = "filtlong_meanq"    
rule quality_filtlong_meanq:    
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq",
		faidx=rules.faidx_self.output.faidx
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log" 
	threads: 8
	resources:
		mem="16GiB",
		runtime=f"{60 * REPEAT}m"
	conda:
		ENVS / "filtlong.yaml"
	params:
		mean_qual="--length_weight 0 --window_q_weight 0",
		target_base=get_target_bases
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.tsv", REPEAT)
	shell:
		"""
		filtlong {params.mean_qual} --target_bases {params.target_base} {input.reads} > {output.reads} 2> {log}
		"""

# 4. Filtlong Length Only (weight 0 for mean_q and window_q)
tool = "filtlong_len"    
rule quality_filtlong_len:    
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq",
		faidx=rules.faidx_self.output.faidx
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log" 
	threads: 8
	resources:
		mem="16GiB",
		runtime=f"{60 * REPEAT}m"
	conda:
		ENVS / "filtlong.yaml"
	params:
		read_length="--mean_q_weight 0 --window_q_weight 0",
		target_base=get_target_bases
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.tsv", REPEAT)
	shell:
		"""
		filtlong {params.read_length} --target_bases {params.target_base} {input.reads} > {output.reads} 2> {log}
		"""

# 5. Filtlong Default (length + quality weighting)
tool = "filtlong_default"    
rule quality_filtlong_default:    
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq",
		faidx=rules.faidx_self.output.faidx
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log" 
	threads: 8
	resources:
		mem="16GiB",
		runtime=f"{60 * REPEAT}m"
	conda:
		ENVS / "filtlong.yaml"
	params:
		default="--length_weight 1 --mean_q_weight 1 --window_q_weight 1",
		target_base=get_target_bases
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.tsv", REPEAT)
	shell:
		"""
		filtlong {params.default} --target_bases {params.target_base} {input.reads} > {output.reads} 2> {log}
		"""

# 6. Chopper Trim-by-Quality 100 bp
tool = "chopper_trim100"
rule quality_chopper_trim100:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq"
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="16GiB",
		runtime=f"{25 * REPEAT}m"
	conda:
		ENVS / "chopper.yaml"
	params:
		approach="--trim-approach trim-by-quality",
		quality_threshold=get_quality_threshold,
		length_threshold=100
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.tsv", REPEAT)
	shell:
		"""
		chopper {params.approach} --cutoff {params.quality_threshold} --minlength {params.length_threshold} -i {input.reads} > {output.reads} 2> {log}
		"""

# Chopper Trim-by-Quality 1000 bp
tool = "chopper_trim1000"
rule quality_chopper_trim1000:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq"
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="16GiB",
		runtime=f"{25 * REPEAT}m"
	conda:
		ENVS / "chopper.yaml"
	params:
		approach="--trim-approach trim-by-quality",
		quality_threshold=get_quality_threshold,
		length_threshold=1000
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.tsv", REPEAT)
	shell:
		"""
		chopper {params.approach} --cutoff {params.quality_threshold} --minlength {params.length_threshold} -i {input.reads} > {output.reads} 2> {log}
		"""

# 7. Chopper Best-Read-Segment Extract 100 bp
tool = "chopper_extract100"
rule quality_chopper_extract100:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq"
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="16GiB",
		runtime=f"{25 * REPEAT}m"
	conda:
		ENVS / "chopper.yaml"
	params:
		approach="--trim-approach best-read-segment",
		quality_threshold=get_quality_threshold,
		length_threshold=100
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.tsv", REPEAT)
	shell:
		"""
		chopper {params.approach} --cutoff {params.quality_threshold} --minlength {params.length_threshold} -i {input.reads} > {output.reads} 2> {log}
		"""

# Chopper Best-Read-Segment Extract 1000 bp
tool = "chopper_extract1000"
rule quality_chopper_extract1000:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq"
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log"
	threads: 8
	resources:
		mem="16GiB",
		runtime=f"{25 * REPEAT}m"
	conda:
		ENVS / "chopper.yaml"
	params:
		approach="--trim-approach best-read-segment",
		quality_threshold=get_quality_threshold,
		length_threshold=1000
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq")
	benchmark:
		repeat(BENCHMARK / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.tsv", REPEAT)
	shell:
		"""
		chopper {params.approach} --cutoff {params.quality_threshold} --minlength {params.length_threshold} -i {input.reads} > {output.reads} 2> {log}
		"""

# 8. Unprocessed Control (no quality filtering)
tool = "unprocessed"
rule quality_notrim:
	input:
		reads=RESULTS / "QC/trimming/{trimmer}/{model}/{sample}.{trimmer}.fastq"
	log:
		LOGS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.log"
	resources:
		mem="8GiB",
		runtime="10m"
	output:
		reads=temp(RESULTS / f"QC/quality/{tool}-{{trimmer}}/{{model}}/{{sample}}.{tool}-{{trimmer}}.fastq")
	shell:
		"""
		cp {input.reads} {output.reads}
		"""
