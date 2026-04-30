process ALLOCOUNT {
	label 'allopipe'
	publishDir "${output_dir}", mode: 'copy', overwrite: false

	input:
	path donor_input
	path recipient_input
	val  run_name
	val  orientation
	val  imputation
	val  optional_args
	val  output_dir

	output:
	path "runs",  emit: results_dir 
	
	script:
	"""
	allopipe_src_dir=${projectDir}/src
	python \${allopipe_src_dir}/ams_pipeline.py \
		-n ${run_name} \
		${donor_input} \
		${recipient_input} \
		${orientation} \
		${imputation} \
		${optional_args} \
		-o ./
	"""
}
