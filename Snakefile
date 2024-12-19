##rules to add:
#adapter_trimming
#fastqc (all non-fastqc-ed files
#demultiplexing

rule bwa_map:
    input:
        "genome/GCA_003113815.1_ASM311381v1_genomic.fna",
        "samples/{sample}.fq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
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
		bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
	output:
		"output_refmap/populations.snps.vcf"
	shell:
		'ref_map.pl --samples sorted_reads --popmap popmap.txt  -o output_refmap -X '"populations:--vcf"''


rule filtering:
	input:
		"output_refmap/populations.snps.vcf"
	output:
		"filtered.recode.vcf"
	shell:
		"vcftools --vcf {input} --max-missing 0.8 --recode --out filtered"
