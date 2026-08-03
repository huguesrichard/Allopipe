process ALLO_COUNT {
	label 'allopipe'
	tag { params.mode == 'cohort' ? pair_id.replaceFirst(/^P/, 'PAIR ') : "PAIR: ${donor_input.simpleName} | ${recipient_input.simpleName}" }
	publishDir "${output_dir}", mode: 'copy', overwrite: true, saveAs: { filename ->
		filename.startsWith('publish/') ? filename.substring('publish/'.length()) : null
	}

	input:
	tuple val(pair_id), path(donor_input), path(recipient_input)
	val   run_name
	val   orientation
	val   imputation
	val   frameshift
	val   allo_count_opts
	val   output_dir
	val   nextflow_command_base64

	output:
	tuple val(pair_id),	val(run_name), path("runs/${run_name}"), emit: results_dir
	path "publish/**"
	
	script:
	"""
	allopipe_src_dir=${projectDir}/src
	export ALLOPIPE_PUBLISHED_OUTPUT_DIR=${output_dir}
	export ALLOPIPE_NEXTFLOW_COMMAND_BASE64=${nextflow_command_base64}
	
	python \${allopipe_src_dir}/ams_pipeline.py \
		-n ${run_name} \
		${donor_input} \
		${recipient_input} \
		${orientation} \
		${imputation} \
		${frameshift ? '--frameshift' : ''} \
		--pair ${pair_id} \
		${allo_count_opts} \
		-o ./

	# Hard-link mirror for file-by-file publishing without replacing the existing run tree
	mkdir -p publish
	cp -al runs publish/
	"""
}
