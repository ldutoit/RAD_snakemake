configfile: "config.yaml"  # Replace with your actual config file name
import os
SAMPLES = []
if os.path.exists("popmap.txt"):
    with open("popmap.txt") as f: #only the ones I want
        for line in f:
            SAMPLES.append(line.split()[0])
else:
	with open("barcodes.txt") as f: #do all the samples
		for line in f:
			SAMPLES.append(line.split()[2])


rule make_popmap:
	input:
		"barcodes.txt"
	output:
		"popmap.txt"
	shell:
		"""
        awk '{{print $3 "\\tpop"}}' {input} > {output}
		"""
if config["mode"] == "denovo":
	rule adapter_trimming:
		input:
			forward_reads=config["raw_fastq"]["forward"],
			reverse_reads=config["raw_fastq"]["reverse"]
		output:
			"trimmed/trimmed_R1_001.fastq",
			"trimmed/trimmed_R2_001.fastq"
		threads: workflow.cores
		shell:
			"""
			mkdir -p trimmed
			cutadapt -j {threads} -a {config[cutadapt][adapter]} -A {config[cutadapt][adapter]} \
			-q {config[cutadapt][minimum_phred]} -o trimmed/trimmed_R1_001.fastq \
			--minimum-length {config[cutadapt][length]}:{config[cutadapt][length]} \
			--length {config[cutadapt][length]} \
			-p trimmed/trimmed_R2_001.fastq {input.forward_reads} {input.reverse_reads}
			fastqc trimmed/trimmed*fastq
			"""
elif config["mode"] == "refmap":
	rule adapter_trimming:
		input:
			forward_reads=config["raw_fastq"]["forward"],
			reverse_reads=config["raw_fastq"]["reverse"]
		output:
			"trimmed/trimmed_R1_001.fastq",
			"trimmed/trimmed_R2_001.fastq"
		threads: workflow.cores
		shell:
			"""
			mkdir -p trimmed
			cutadapt -j {threads} -a {config[cutadapt][adapter]} -A {config[cutadapt][adapter]} \
			-q {config[cutadapt][minimum_phred]} -o trimmed/trimmed_R1_001.fastq \
			--minimum-length {config[cutadapt][length]}:{config[cutadapt][length]} \
			-p trimmed/trimmed_R2_001.fastq {input.forward_reads} {input.reverse_reads}
			fastqc trimmed/trimmed*fastq
			"""
rule demultiplex:
    input:
        r1="trimmed/trimmed_R1_001.fastq",
        r2="trimmed/trimmed_R2_001.fastq",
        barcodes="barcodes.txt"
    output:
        fqgz1=expand("samples/{sample}.1.fq.gz", sample=SAMPLES),
        fqgz2=expand("samples/{sample}.2.fq.gz", sample=SAMPLES)
    threads: workflow.cores  # Adjust based on your requirements
    shell:
        """
        mkdir -p raw
		cd raw
		ln -sf ../{input.r1} .
        ln -sf ../{input.r2} .
		cd ..
        # Running process_radtags with provided parameters
        process_radtags -P -p raw/ -o ./samples/ -b {input.barcodes} -e pstI -r -c -q --inline-inline
        """
rule index_genome:
	input:
		genome=config["genome"]["ref"]
	output:
		"{ref}.amb"
	shell:
		"bwa index {genome}"

rule bwa_map:
    input:
        genome=config["genome"]["ref"],
        read1="samples/{sample}.1.fq.gz",
        read2="samples/{sample}.2.fq.gz"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input.genome} {input.read1} {input.read2} | samtools view -Sb - > {output}"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

if config["mode"] == "denovo":
	rule variant_calling:
		input:
			fqgz1=expand("samples/{sample}.1.fq.gz", sample=SAMPLES),
			fqgz2=expand("samples/{sample}.2.fq.gz", sample=SAMPLES),
			popmap="popmap.txt"
		output:
			"output_SNPcalling/populations.snps.vcf"
		shell:
			'denovo_map.pl --paired --samples samples --popmap popmap.txt -o output_SNPcalling -X '"populations:--vcf"''
elif config["mode"] == "refmap":
	rule variant_calling:
		input:
			bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
			popmap="popmap.txt"
		output:
			"output_SNPcalling/populations.snps.vcf"
		threads: workflow.cores		
		shell:
			'ref_map.pl --samples sorted_reads --popmap popmap.txt -o output_SNPcalling -X '"populations:--vcf"''


rule filtering:
	input:
		vcf="output_SNPcalling/populations.snps.vcf",
	output:
		"filtered.recode.vcf",
		"filtered.imiss"
	threads: workflow.cores
	shell:
		"""
		vcftools --vcf {input.vcf} {config[vcf_filtering][parameters]} --recode --out filtered
		vcftools --vcf filtered.recode.vcf --missing-indv --out filtered
		"""

