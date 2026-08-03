process ALLO_AFFINITY {
	label 'allopipe'
	tag { params.mode == 'cohort' ? pair_id.replaceFirst(/^P/, 'PAIR ') : "PAIR: ${donor_input.simpleName} | ${recipient_input.simpleName}" }
	publishDir "${output_dir}", mode: 'copy', overwrite: true, saveAs: { filename ->
		filename.startsWith('publish/') ? filename.substring('publish/'.length()) : null
	}

	input:
	tuple val(pair_id), path(donor_input), path(recipient_input), val(hla_typing), val(run_name), path(allo_count_results)
	val   ensembl_path
	val   allo_affinity_opts
	val   output_dir

	output:
	tuple val(pair_id), val(run_name), path("runs/${run_name}"), emit: results_dir
	path "publish/**"
	
	script:
	"""
	mkdir -p runs/${run_name}
	cp -a ${allo_count_results}/. runs/${run_name}/

	allopipe_src_dir=${projectDir}/src
	python \${allopipe_src_dir}/aams_pipeline.py \
		-n ${run_name} \
		-p ${pair_id} \
		-d ${projectDir}/${ensembl_path} \
		-a ${hla_typing} \
		${allo_affinity_opts} \
		-o ./

	# Hard-link mirror for file-by-file publishing without replacing the existing run tree
	mkdir -p publish
	cp -al runs publish/
	"""
}
