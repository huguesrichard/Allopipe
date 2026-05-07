process ALLO_COUNT {
	label 'allopipe'
	tag "$pair_id"
	publishDir "${output_dir}", mode: 'copy', overwrite: false, enabled: params.mode == 'pair'

	input:
	tuple val(pair_id), path(donor_input), path(recipient_input)
	val  run_name
	val  orientation
	val  imputation
	val  allo_count_opts
	val  output_dir

	output:
	tuple val(pair_id), val(run_name), path("runs/${run_name}"), emit: results_dir
	
	script:
	"""
	allopipe_src_dir=${projectDir}/src
	python \${allopipe_src_dir}/ams_pipeline.py \
		-n ${run_name} \
		${donor_input} \
		${recipient_input} \
		${orientation} \
		${imputation} \
		--pair ${pair_id} \
		${allo_count_opts} \
		-o ./
	"""
}
