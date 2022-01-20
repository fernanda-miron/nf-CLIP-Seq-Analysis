#!/usr/bin/env nextflow

/*================================================================
The Aguilar Lab presents...

The CLIP-Seq analysis Pipeline

- A tool for CLIP-Seq data analysis from pre-processing to
peak-calling

==================================================================
Version: 0.1
Project repository:
==================================================================
Authors:

- Bioinformatics Design
 Diana Rogel (diana.rogel@mpi-bn.mpg.de)
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)
 Fernanda Miron-Toruno (fernandamiront@gmail.com)

- Bioinformatics Development
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)
 Fernanda Miron-Toruno (fernandamiront@gmail.com)

- Nextflow Port
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)
 Fernanda Miron-Toruno (fernandamiront@gmail.com)

=============================
Pipeline Processes In Brief:

Pre-processing:
- fastq file alignment
- sam to bam transformation

Core-processing:
- peak calling

Pos-processing:
- bed annotation

Analysis:
================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
	The CLIP-Seq Pipeline
  v${version}
  ==========================================

	Usage:

	nextflow run ${pipeline_name}.nf --fastq_file <path to input 1> [--output_dir path to results ]

	  --fastq_file	<- fastq.gz fasta file;

	  --output_dir     <- directory where results, intermediate and log files will bestored;
	      default: same dir where --query_fasta resides

	  -resume	   <- Use cached results if the executed project has been run before;
	      default: not activated
	      This native NF option checks if anything has changed from a previous pipeline execution.
	      Then, it resumes the run from the last successful stage.
	      i.e. If for some reason your previous run got interrupted,
	      running the -resume option will take it from the last successful pipeline stage
	      instead of starting over
	      Read more here: https://www.nextflow.io/docs/latest/getstarted.html#getstart-resume
	  --help           <- Shows Pipeline Information
	  --version        <- Show version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "nf-CLIP-Seq"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.fastq_file = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "${pipeline_name} v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at MAY 2021
*/
nextflow_required_version = '20.01.0'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*//////////////////////////////
  INPUT PARAMETER VALIDATION BLOCK
*/

/* Check if the input directory is provided
    if it was not provided, it keeps the 'false' value assigned in the parameter initiation block above
    and this test fails
*/
if ( !params.fastq_file) {
  log.error " Please provide the --fastq_file \n\n" +
  " For more information, execute: nextflow run compare-miRNA-pairs.nf --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --input_dir
*/
params.output_dir = file(params.fastq_file).getParent() //!! maybe creates bug, should check

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable pipeline_name defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*
Useful functions definition
*/

/*//////////////////////////////
  LOG RUN INFORMATION
*/
log.info"""
==========================================
The CLIP-Seq Analysis Pipeline
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['Input fastq.gz file']			= params.fastq_file
pipelinesummary['Results Dir']		= results_dir
pipelinesummary['Intermediate Dir']		= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

/*//////////////////////////////
  PIPELINE START
*/

/*
	READ INPUTS
*/

/* Load fastq file into channel */
Channel
	.fromPath( "${params.fastq_file}" )
	// .view()
	.set{ clip_input}

/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mk_modules/mk_alignment/*")
	.toList()
	.set{ mkfiles_pre1 }

/* Load reference_index*/
Channel
	.fromPath("${workflow.projectDir}/test/data/reference-index")
	.set{ reference_index}

/* Writting process */
process mk_alignment {

	publishDir "${results_dir}/mk_alignment/",mode:"copy"

	input:
	file fastq from clip_input
	path index from reference_index
	file mk_files from mkfiles_pre1

	output:
	file "*.sam" into results_mk_alignment

	"""
	bash runmk.sh
	"""

}

/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mk_modules/mk_sam_to_bam/*")
	.toList()
	.set{ mkfiles_pre2 }

/* Writting process */
process mk_sam_to_bam {

		publishDir "${results_dir}/mk_sam_to_bam/",mode:"copy"

		input:
		file sam from results_mk_alignment
		file mk_files from mkfiles_pre2

		output:
		file "*.bam" into results_mk_sam_to_bam_bam
		file "*.bai" into results_mk_sam_to_bam_bai

		"""
		bash runmk.sh
		"""

	}

/* Read mkfile module files */
	Channel
		.fromPath("${workflow.projectDir}/mk_modules/mk_peak_calling/*")
		.toList()
		.set{ mkfiles_core1 }

/* Writting process */
process mk_peak_calling {

				publishDir "${results_dir}/mk_peak_calling/",mode:"copy"

				input:
				file bam from results_mk_sam_to_bam_bam
				file bai from results_mk_sam_to_bam_bai
				file mk_files from mkfiles_core1

				output:
				file "*.bed" into results_mk_peak_calling

				"""
				bash runmk.sh
				"""

			}
