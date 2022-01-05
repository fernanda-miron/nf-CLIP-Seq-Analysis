#!/usr/bin/env nextflow

/*Modulos de prueba para  nextflow */
results_dir = "./test/results"
intermediates_dir = "./test/results/intermediates"

process making_index {

publishDir "${results_dir}/reference-index/", mode:"copy"

	input:
	file reference_genome

	output:
	path "*.bt2"

	"""
	bowtie2-build ${reference_genome} reference-index
	"""
}

process making_alignment {

publishDir "${results_dir}/alignment-results/", mode:"copy"

	input:
	path p1
	file gz_file

	output:
	path "*.bam"

	"""
	mkdir reference-index
	mv ${p1} reference-index
	bowtie2 -x reference-index/reference-index -U ${gz_file} | samtools view -bS -> out.bam
	"""
}
