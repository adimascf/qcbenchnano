from pathlib import Path

def get_original_fastqs(wildcards, trim_status):
	# trim_status expects "trim" or "notrim"
	reads_dir = Path(pep.get_sample(wildcards.sample)["reads_dir"])
	return (reads_dir / f"{trim_status}/{wildcards.model}/{wildcards.sample}.fastq")

def get_reference_genome(wildcards):
	ref_dir = Path(pep.get_sample(wildcards.sample)["reference_path"])
	return (ref_dir / f"{wildcards.sample}.fna")

def get_mutreference_genome(wildcards):
	ref_dir = Path(pep.get_sample(wildcards.sample)["mutreference_path"])
	return (ref_dir / f"{wildcards.sample}.mutref.fna")

def get_truth_vcf(wildcards):
	ref_dir = Path(pep.get_sample(wildcards.sample)["truth_vcf_path"])
	return (ref_dir / f"{wildcards.sample}.truth.vcf.gz")

def get_sequencing_kits(wildcards):
	return pep.get_sample(wildcards.sample)["sequencing_kits"]

def get_quality_threshold(wildcards):
	if wildcards.model == "hac":
		return 10
	elif wildcards.model == "sup":
		return 15
	else:
		return 15 # default

def get_target_bases(wildcards, input):
	# sum the lengths of all contigs (column 2 in the .fai file)
	with open(input.faidx, "r") as lines:
		genome_size = sum(int(line.split("\t")[1]) for line in lines)
	return genome_size * 100
