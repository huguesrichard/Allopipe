process ALLOCOUNT {
	label 'allopipe'
	//stageInMode 'copy'

	input:
	path donor_input
	path recipient_input
	val  run_name
	val  orientation
	val  imputation
	val  optional_args

	//output:
	
	script:
	"""
	allopipe_src_dir=${projectDir}/src
	python \${allopipe_src_dir}/ams_pipeline.py -n ${run_name} ${donor_input} ${recipient_input} ${orientation} ${imputation} ${optional_args}
	"""
}