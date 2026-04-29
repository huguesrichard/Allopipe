process ALLOAFFINITY {
	label 'allopipe'
	//stageInMode 'copy'

	input:
	run_name
	ensembl_path
	hla_typing
	optional_args

	//output:
	
	script:
	"""
	allopipe_src_dir=${projectDir}/src
	python \${allopipe_src_dir}/aams_pipeline.py -n ${run_name} -d ${ensembl_path} -a ${hla_typing} ${optional_args} --dry_run
	"""
}