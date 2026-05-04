process ALLOAFFINITY {
	label 'allopipe'
	publishDir "${output_dir}", mode: 'copy', overwrite: true

	input:
	path allocount_results
	val  run_name
	val  ensembl_path
	val  hla_typing
	val  alloaffinity_opts
	val  output_dir

	output:
    path "runs/${run_name}",  emit: affinity_results
	
	script:
	"""
	allopipe_src_dir=${projectDir}/src
	python \${allopipe_src_dir}/aams_pipeline.py \
		-n ${run_name} \
		-d ${projectDir}/${ensembl_path} \
		-a ${hla_typing} \
		${alloaffinity_opts} \
		-o ./ \
		--dry_run
	"""
}
