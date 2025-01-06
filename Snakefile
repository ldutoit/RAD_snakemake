configfile: "config.yaml"
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

rule adapter_trimming:
    input:
        forward_reads=config["raw_fastq"]["forward"],
        reverse_reads=config["raw_fastq"]["reverse"]
    output:
        "trimmed/trimmed_R1_001.fastq",
        "trimmed/trimmed_R2_001.fastq"
    shell:
        """
        mkdir -p trimmed
        cutadapt -a {config[cutadapt][adapter]} -A {config[cutadapt][adapter]} \
        -q {config[cutadapt][minimum_phred]} -o trimmed/trimmed_R1_001.fastq \
        --minimum-length {config[cutadapt][minimum_length]}:{config[cutadapt][minimum_length]} \
        -p trimmed/trimmed_R2_001.fastq {input.forward_reads} {input.reverse_reads}
        fastqc trimmed/trimmed*fastq
        """

rule demultiplex:
	input:
		"trimmed/trimmed_R1_001.fastq",
		"trimmed/trimmed_R2_001.fastq"
	output:
		fqgz1=expand("samples/{sample}.1.fq.gz", sample=SAMPLES),
		fqgz2=expand("samples/{sample}.2.fq.gz", sample=SAMPLES),
	shell:
        	"""
	        mkdir -p raw;
	        ln -sf $PWD/trimmed/*fastq $PWD/raw/;
	        process_radtags -P -p raw/ -o ./samples/ -b barcodes.txt -e pstI -r -c -q --inline-inline
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

rule refmap_variant_calling:
	input:
		bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
		popmap="popmap.txt"
	output:
		"output_refmap/populations.snps.vcf"
	shell:
		'ref_map.pl --samples sorted_reads --popmap popmap.txt -o output_refmap -X '"populations:--vcf"''


rule filtering:
	input:
		vcf="output_refmap/populations.snps.vcf",
	output:
		"filtered.recode.vcf"
	shell:
		"vcftools --vcf {input.vcf} {config[vcf_filtering][parameters]} --recode --out filtered"
