#!/usr/bin/env nextflow

/*Modulos de prueba para  nextflow */
results_dir = "./test/results"
intermediates_dir = "./test/results/intermediates"

process making_index {

	input:
	file reference_genome

	output:
	path "*.bt2"

	"""
	bowtie2-build ${reference_genome} reference-index
	"""
}

process making_alignment {

	input:
	path p1
	file gz_file

	output:
	path "*.sam"

	"""
	mkdir reference-index
	mv ${p1} reference-index
	bowtie2 -x reference-index/reference-index -U ${gz_file} -S out.sam
	"""
}

process sam_to_bam {

	input:
	file p2

	output:
	path "*.bam"

	"""
	samtools view -S -b ${p2} > out.bam
	"""
}

process sorting_bam {

	input:
	file p3

	output:
	path "*.bam"

	"""
	samtools sort ${p3} -o sample.sorted.bam
	"""
}

process peak_calling {
	publishDir "${results_dir}/peak-call/", mode:"copy"

	input:
	file p4

	output:
	path "*.bed"

	"""
	source activate clipper3
	clipper -b ${p4} -o out.bed -s GRCh38
	"""
}

process peak_calling {
	publishDir "${results_dir}/peak-call/", mode:"copy"

	input:
	file p4

	output:
	path "*.bed"

	"""
	source activate clipper3
	clipper -b ${p4} -o out.bed -s GRCh38
	"""
}
